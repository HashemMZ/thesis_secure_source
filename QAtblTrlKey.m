function [code,codeh,ST]=QAtblTrlKey(seq,trellis,N,key)
%encoding using trellis of code
%key=[0  0  0  28  0    0   0 ]
%ST  =[1 22 3  14  X    X    X]
%ST =[              28  2   11  29]
len=length(seq);
code=[];codeh=[];
%make lookup tables for N
% seq(find(seq==1) )=0;
% seq(find(seq>1) )=1;
ST=key(1);
t=1;
while t<=len
    sym=seq(t);
    ok=0;
    while ~ok
        for i=1: trellis(ST(end)).outNo
            if isequal(sym,trellis(ST(end)).in(i).code)%find proper output branch of trellis                
                R=trellis(ST(end)).out(i).code;
                Rh=trellis(ST(end)).Huffout(i).code;
                ST=[ST trellis(ST(end)).outstate(i)];
                if key(length(ST))
                    ST(end)=key(length(ST));
                end
                codeh=[codeh Rh];
                code=[code R];
                ok=1;
                break
            end
        end
        clc
        t=t+1
        pause(0.00001)
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
    Rh=trellis(ST(end)).Huffout(1).code;
    ST=[ST trellis(ST(end)).outstate(1)];
    code=[code R];
    codeh=[codeh Rh];
    rembits=rembits-length(R);
end
    
