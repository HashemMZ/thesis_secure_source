function seq=QADtblTrlKeyH(codeh,trellis,len,key)
%decoding using huffman output trellis of code
ST=key(1);seq=[];
lenC=length(codeh);
t=1;
while length(seq)<len %until decoding last bit
    last=1;
    for i=1: trellis(ST(end)).outNo % get maximum length of Huffman output in trellis for current state,ex: if outputs are 11-01-00-101-100 then max length is 3
        last=max(last,length(trellis(ST(end)).Huffout(i).code));
    end
    if t+last-1>length(codeh)
        last=length(codeh)-t+1;
    end
    out=codeh(t:t+last-1);
    ok=0;
    while ~ok
        for i=1: trellis(ST(end)).outNo
            if isequal(out,trellis(ST(end)).Huffout(i).code)%find proper output branch of trellis
                dec=trellis(ST(end)).in(i).code;
                ST=[ST trellis(ST(end)).outstate(i)];
                if key(length(ST))
                    ST(end)=key(length(ST));
                end
                seq=[seq dec];
                ok=1;
                t=t+length(out)
                break
            end
        end

        if ~ok %if proper output branch needs less output bits
            out(end)=[];
        end
    end
end
seq=seq(1:len);