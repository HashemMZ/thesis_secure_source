function crc=CRC(msg,gen)
% compute CRC for message with generator polynomial gen
% uses systematic form (i.e msg bits and check bits are appended )
%in inserting msg and gen left most bit is the coefficient of highest degree ex: X^3 + X + 1 =[1 0 1 1]
% if no gen is given we use X.25 CRC generator
%Hashem Moradmand 2008/8/11 ,Persian 1387/5/21
if nargin==1%X.25 CRC=X^16+X^12+X^5+1
    gen=[1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1];
end
k=length(msg);
r=length(gen);
n=k+r-1;% r=n-k check bits
% multiply msg with X^(r-1)
msg=[msg zeros(1,r-1)];
%Division by gen
remainder=msg(1:r);
for pnt=r:n      
    if remainder(1)      
            remainder=bitxor(remainder,gen);      
    end    
    if pnt<n
        %shift rem%next bit
        remainder=[remainder(2:end) msg(pnt+1)];        
    end
end
    
%check bits from remainder
 crc=remainder(2:end);



