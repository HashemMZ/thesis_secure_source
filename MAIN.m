%MAIN combined (SJSC)
%--------------------------------------------------------------------------
clear all
format long
while 1
    clc
    reply = input('Do you want continue last simulation? y/n [y]: ', 's');
    if isempty(reply) || reply=='y'
        param=input('param? ');
        snr=input('snr? ');
        eval( ['load  F:\VCRC\',int2str(param),'_',num2str(snr),'\loader']);
        disp(['now siml=',num2str(siml)])
        stp=input('stop? ');        
        strt=siml+1;        
        break
    elseif reply == 'n';
        strt=input('start? ');
        stp=input('stop? ');
        param=input('param? ');
        snr=input('snr? ');
        every=input('every? ');
        cBlk=0;cBit=0;eBlk=0;eBit=0;PER=0;BER=0;SER=0;eSym=0;cSym=0;rate=0;
        ccBlk=0;ccBit=0;ceBlk=0;ceBit=0;cPER=0;cBER=0;cSER=0;ceSym=0;ccSym=0;
        break
    end
end

% eval( ['load F:\VCRC',int2str(param),'\est',int2str(param)']);
eval( ['load  F:\VCRC\',int2str(param),'_',num2str(snr),'\enctbl']);
trelHuffout
% t=zeros(1,stp-strt+1);
%88888888888888888888888888888888888888888888 start simulation
for siml=strt:stp

    disp(['simulation: ',num2str(siml),' every=',num2str(every),' param=',num2str(param),' SNR =',num2str(snr),])
    disp(['BER =',num2str(BER),' PER =',num2str(PER),' SER =',num2str(SER),'//', ' cBER =',num2str(cBER),' cPER =',num2str(cPER),' cSER =',num2str(cSER)])
%     tic
    R=[];len=256;
    seq=randsrc(1,len,[1 3;counts(1,[1 3])./sum(counts(1,[1 3]))]);%seq=randsrc(1,blklen,[[1 2]; prbs]) ;
    crc=CRC((seq-1)/2);
    seqcrc=[1+2*crc seq];
    lencrc=length(seqcrc);
    
    place = randsrc(1,lencrc,[ 0 1 ; 0.5 0.5]);%!
    place(1)=1;
    stateperm=randsrc(1,lencrc,[1:sNo]);
    key=stateperm.*place;   
    [codes,R,ST]=QAtblTrlKey(seqcrc,trellis,N,key);      
%     [R,ST]=QAtblTrl(seqcrc,trellis,N); % <<<<<<<<<<<<<<<<<<<<KEY   , maybe a key should be generated each time 
    code= awgn(R*2-1,snr,'measured');
    %     rsrvCode=code;
    r_hrd=(sign(code)+1)/2;
    errpoints=find ( r_hrd ~= R );
    lenC=length(code)+1;%Viterbi needs
    rate=rate+(len/(lenC-1));
    %88888888888888888888888888888888888888888%  Viterbi
    L=10;simpVit=0;
    format;
    [highcost,s_info]=LVDC1(R,code,trellis,L,simpVit,key,ST);%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<chng
    format long;
    crcOK=0;
    for pth=1:L
        [dseq bitclk]= sqarithdecoflush(highcost(pth).path,counts,N,Fmax,lencrc,midFS);%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<chng
        if ~ismember(2,dseq)
            if isequal(dseq(1:16),1+2*CRC((dseq(17:end)-1)/2)) && isequal(dseq(17:end),seq)
                if pth==1%if it is not first path it means that crc has detected this point
                    ccBlk=ccBlk+1;%add correct crc
                    ccBit=ccBit+lenC-1;
                    ccSym=cSym+len;
                    cPER=ceBlk/(ceBlk+ccBlk);
                    cSER=ceSym/(ceSym+ccSym);
                    cBER=ceBit/(ceBit+ccBit);
                    %88888888888888888888888
                    cBlk=cBlk+1;%add correct 
                    cBit=cBit+lenC-1;
                    cSym=cSym+len;
                    PER=eBlk/(eBlk+cBlk);
                    SER=eSym/(eSym+cSym);
                    BER=eBit/(eBit+cBit);
                else
                    ccBlk=ccBlk+1;%add correct crc
                    ccBit=ccBit+lenC-1;
                    ccSym=cSym+len;
                    cPER=ceBlk/(ceBlk+ccBlk);
                    cSER=ceSym/(ceSym+ccSym);
                    cBER=ceBit/(ceBit+ccBit);
                    %88888888888888888888888error other
                    [dseq bitclk]= sqarithdecoflush(highcost(1).path,counts,N,Fmax,lencrc,midFS);
                    eSym=eSym+sum((dseq(17:end)~=seq));
                    cSym=cSym+len-sum((dseq(17:end)~=seq));
                    eBit=eBit+sum(xor(highcost(1).path,R));
                    cBit=cBit+lenC-1-sum(xor(highcost(1).path,R));
                    eBlk=eBlk+1;
                    PER=eBlk/(eBlk+cBlk);
                    BER=eBit/(eBit+cBit);
                    SER=eSym/(eSym+cSym);
                end
%                 t(siml) = toc;
                if rem(siml,every)==0
                    DATA(siml/every,:)=[cBlk,cBit,cSym,eBlk,eBit,eSym,PER,BER,SER,ccBlk,ccBit,ccSym,ceBlk,ceBit,ceSym,cPER,cBER,cSER];
                    eval( ['save F:\VCRC\',int2str(param),'_',num2str(snr),'\s', num2str(siml),' DATA']);
                    eval( ['save F:\VCRC\',int2str(param),'_',num2str(snr),'\loader',' siml',' cBlk',' cBit',' eBlk',' eBit',' eSym',' cSym',' PER',' BER',' SER',' ccBlk',' ccBit',' ceBlk',' ceBit',' ceSym',' ccSym',' cPER',' cBER',' cSER',' param',' every',' stp',' DATA',' rate']);
                end
                crcOK=1;
                break;
            end
        end
    end
    if crcOK
        continue;
    end
    
    %888888888888888888888888888888888888888888888888888888 decoded seq || CRC not OK
    [dseq bitclk]= sqarithdecoflush(highcost(1).path,counts,N,Fmax,lencrc,midFS);
    eSym=eSym+sum((dseq(17:end)~=seq));
    cSym=cSym+len-sum((dseq(17:end)~=seq));
    eBit=eBit+sum(xor(highcost(1).path,R));
    cBit=cBit+lenC-1-sum(xor(highcost(1).path,R));
    eBlk=eBlk+1;
    PER=eBlk/(eBlk+cBlk);
    BER=eBit/(eBit+cBit);
    SER=eSym/(eSym+cSym);
    %8888888888888888888888 error crc 
    ceSym=ceSym+sum((dseq(17:end)~=seq));
    ccSym=ccSym+len-sum((dseq(17:end)~=seq));
    ceBit=ceBit+sum(xor(highcost(1).path,R));
    ccBit=ccBit+lenC-1-sum(xor(highcost(1).path,R));
    ceBlk=ceBlk+1;
    cPER=ceBlk/(ceBlk+ccBlk);
    cBER=ceBit/(ceBit+ccBit);
    cSER=ceSym/(ceSym+ccSym);
    eval( ['save  F:\VCRC\',int2str(param),'_',num2str(snr),'\decfail', int2str(siml)]);
    
    %     clc,siml,sum(dseq-seq) ,pause(.000001)
    %     t(siml) = toc;
    if rem(siml,every)==0
        DATA(siml/every,:)=[cBlk,cBit,cSym,eBlk,eBit,eSym,PER,BER,SER,ccBlk,ccBit,ccSym,ceBlk,ceBit,ceSym,cPER,cBER,cSER];
        eval( ['save F:\VCRC\',int2str(param),'_',num2str(snr),'\s', num2str(siml),' DATA']);
        eval( ['save F:\VCRC\',int2str(param),'_',num2str(snr),'\loader',' siml',' cBlk',' cBit',' eBlk',' eBit',' eSym',' cSym',' PER',' BER',' SER',' ccBlk',' ccBit',' ceBlk',' ceBit',' ceSym',' ccSym',' cPER',' cBER',' cSER',' param',' every',' stp',' DATA',' rate']);
    end

end