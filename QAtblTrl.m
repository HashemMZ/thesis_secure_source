function [code,ST]=QAtblTrl(seq,trellis,N)
%encoding using trellis of code
len=length(seq);
code=[];
%make lookup tables for N
% seq(find(seq==1) )=0;
% seq(find(seq>1) )=1;
ST=1;
t=1;
while t<=len
    sym=seq(t);
    ok=0;
    while ~ok
        for i=1: trellis(ST(end)).outNo
            if isequal(sym,trellis(ST(end)).in(i).code)%find proper output branch of trellis                
                R=trellis(ST(end)).out(i).code;
                ST=[ST trellis(ST(end)).outstate(i)];
                code=[code R];
                ok=1;
                break
            end
        end
        t=t+1;
        if ~ok %if proper output branch needs more input bits                       
            if t>len
                sym=[sym 1];%            
            else
                sym=[sym seq(t)];%
            end
        end
    end   
end

% Terminate encoding
rembits=N+trellis(ST(end)).fol;%N+E3_count bits is needed to finish encoding
while rembits>0
    R=trellis(ST(end)).out(1).code;%we encode dummy inputs to flush out N bits
    ST=[ST trellis(ST(end)).outstate(1)];
    code=[code R];
    rembits=rembits-length(R);
end
    
