function [code,state] = sqarithencoflush(symb,counts,state,N,E3_max,midFS,EOB)
%
%    [code,state] = sqarithencoflush(symb,counts,state,N,E3_max,midFS)
%
% Encode a symbol using quasi-arithmetic coding. it is
% State oriented, in Renorm process it limits follow to E3_max,N is full
% interval length, state=[low high follow] in which high is closed end,
% counts is the probablities or counts of symbols, set midFS if the
% FS pattern is (a FS b) and clear it if (a b FS)
dec_low=state(1);dec_up=state(2);E3_count=state(3);
if midFS% check the FS pattern (a FS b) or (a b FS)
    A=1;B=3;
else
    A=1;B=2;
end
    
% Compute the cumulative counts vector from the counts
cum_counts = [0, cumsum(counts)];
total_count = cum_counts(end);
code=[];
HALF=2^N/2;
code_index = 1;
symbol = symb;


%     Compute the new  bound
dec_low_new = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1-1)/total_count );
dec_up = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1)/total_count )-1;
dec_low = dec_low_new;


% Check for E1, E2 or E3 conditions and keep looping as long as they occur.
while( isequal(bitget(dec_low, N), bitget(dec_up, N)) || ...
        (isequal(bitget(dec_low, N-1), 1) && isequal(bitget(dec_up, N-1), 0) ) ),
    if (E3_count>=E3_max) && (dec_up >= HALF &&  dec_up <1.5* HALF && dec_low < HALF && dec_low >= HALF/2)
        % Compute the cumulative counts vector from the counts
        if symbol==B
            dec_low=HALF;
        elseif symbol==A
            dec_up=HALF-1;
        end
    end
    % If it is an E1 or E2 condition,
    if isequal(bitget(dec_low, N), bitget(dec_up, N)),
        b = bitget(dec_low, N);
        code(code_index) = b;
        code_index = code_index + 1;
        dec_low = bitshift(dec_low, 1) + 0;
        dec_up = bitshift(dec_up, 1) + 1;
        if (E3_count > 0),
            code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
            code_index = code_index + E3_count;
            E3_count = 0;
        end
        dec_low = bitset(dec_low, N+1, 0);
        dec_up  = bitset(dec_up, N+1, 0);
    elseif ( (isequal(bitget(dec_low, N-1), 1) && ...
            isequal(bitget(dec_up, N-1), 0) ) ),
        dec_low = bitshift(dec_low, 1) + 0;
        dec_up  = bitshift(dec_up, 1) + 1;
        dec_low = bitset(dec_low, N+1, 0);
        dec_up  = bitset(dec_up, N+1, 0);
        dec_low = bitxor(dec_low, 2^(N-1) );
        dec_up  = bitxor(dec_up, 2^(N-1) );
        E3_count = E3_count+1;
    end
end


if(EOB)
    % Terminate encoding
    bin_low = de2bi(dec_low, N, 'left-msb');
    if E3_count==0,
        % Just transmit the final value of the lower bound bin_low
        code(code_index:code_index + N - 1) = bin_low;
        code_index = code_index + N;
    else
        % Transmit the MSB of bin_low.
        b = bin_low(1);
        code(code_index) = b;
        code_index = code_index + 1;

        % Then transmit complement of b (MSB of bin_low), E3_count times.
        code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
        code_index = code_index + E3_count;

        % Then transmit the remaining bits of bin_low
        code(code_index:code_index+N-2) = bin_low(2:N);
        code_index = code_index + N - 1;
    end
end
state=[dec_low,dec_up,E3_count];
% Output only the filled values
code = code(1:code_index-1);