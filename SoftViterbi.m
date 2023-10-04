clear all
load start2
% load enctbl
% R=[];len=256;snr=5;midFS=1;
% seq=randsrc(1,len,[1 3;counts(1,[1 3])./sum(counts(1,[1 3]))]);%seq=randsrc(1,blklen,[[1 2]; prbs]) ;
% crc=CRC((seq-1)/2);crc(find(crc==1))=1+2*midFS;crc(find(crc==0))=midFS;
% seqcrc=[crc seq];
% lencrc=length(seqcrc);
% 
% place = randsrc(1,lencrc,[ 0 1 ; 0.9 0.1]);%!
% place(1)=1;
% stateperm=randsrc(1,lencrc,[1:sNo]);
% key=stateperm.*place;   
% % key=zeros(size(key));key(1)=1;
% [codes,R,ST]=QAtblTrlKey(seqcrc,trellis,N,key);
% dseq=QADtblTrlKeyH(R,trellis,lencrc,key)
% %     [R,ST]=QAtblTrl(seqcrc,trellis,N); % <<<<<<<<<<<<<<<<<<<<KEY   , maybe a key should be generated each time
% code= awgn(R*2-1,snr,'measured');
% %     rsrvCode=code;
% r_hrd=(sign(code)+1)/2;
% errpoints=find ( r_hrd ~= R );
% lenC=length(code)+1;%Viterbi needs
%88888888888888888888888888888888888888888%  Viterbi
L=10;simpVit=0;
clear s_info
lenC=length(code)+1;%code frame bit length
sNo=length(trellis);
% %filling the state information struct
s_info= struct('cost',0,'path',[],'NoPath',0,'intime',0,'hasin',0); %'instates',[],
s_info.path=struct('dec',[],'ST',[]);
s_info=repmat(s_info,sNo,lenC);

%initial state in  trellis is t=1 NOT t=0
t=1;firstst=key(1);
% s_info(firstst,t).instates=0;
s_info(firstst,t).path.dec=[];
s_info(firstst,t).path.ST=firstst;%<<<<<<<<<<<<<<<<
s_info(firstst,t).NoPath=1;
s_info(firstst,t).hasin=1;
disp(['viterbi L=',num2str(L)])
%viterbi
for t=1:lenC
    disp(['  t= ',int2str(t)])
    chng=[];
    for s_i=1:sNo%                                       k--------i------------j           for i (for j (for k)))
        %         disp(['  t= ',int2str(t),'  st= ',int2str(s_i)])
        %88888888888888888888888888888888888 sorting inputs to s_i and restrict them to L ones
        %         tmpstate=s_info(s_i,t).instates;
        tmpcost=s_info(s_i,t).cost;
        tmptime=s_info(s_i,t).intime;
        for q=1:s_info(s_i,t).NoPath
            tmpdecoded.path(q).dec=s_info(s_i,t).path(q).dec;
            tmpdecoded.path(q).ST=s_info(s_i,t).path(q).ST;
        end
        [sortedC,map]=sort(tmpcost,'descend');%sorting the costs
        for in=1:min(L, s_info(s_i,t).NoPath)
            %             s_info(s_i,t).instates(in)=tmpstate(map(in));
            s_info(s_i,t).intime(in)=tmptime(map(in));
            s_info(s_i,t).cost(in)=sortedC(in);
            s_info(s_i,t).path(in).dec=tmpdecoded.path(map(in)).dec;
            s_info(s_i,t).path(in).ST=tmpdecoded.path(map(in)).ST;%<<<<<<<<<<<<<<<
%             s_info(s_i,t).NoPath=in;
        end
        if L< s_info(s_i,t).NoPath% constrain input numbers to L ones
            %             s_info(s_i,t).instates(L+1:end)=[];
            s_info(s_i,t).intime(L+1:end)=[];
            s_info(s_i,t).cost(L+1:end)=[];
            s_info(s_i,t).path(L+1:end)=[];
            s_info(s_i,t).NoPath=L;
        end
        % no need to check hasin flag?<<<<<<<<<<<<<<<<<<<<<<<<<
        if s_info(s_i,t).hasin~=1%if this state is not passed before, thus its outputs are not important
            continue
        end
        %888888888888888888888888888888888888888888888 output states
        for j=1:trellis(s_i).outNo
            codlenH=length(trellis(s_i).Huffout(j).code);%Huffman code length
            s_j=trellis(s_i).outstate(j);tj=t+codlenH;codej= trellis(s_i).Huffout(j).code;% select output states s_j from state s_i
            %             [tj,s_j]
            if tj>lenC || (tj==lenC && (s_j~=trellis(s_i).outstate(1) || tj~=t+codlenH))%ignore final transitions that will go further than the whole frame OR don't obey the tail bitting
                continue
            end

            %88888888888888888888888888888888888 Input states
            inNo=s_info(s_i,t).NoPath;
            for k=1:inNo-simpVit*(inNo-1)
                nexstNo=length(s_info(s_i,t).path(k).ST)+1;
                jumpst=key(nexstNo);
                if jumpst  %if key imposes a jump (not zero )
                    ind=s_info(jumpst,tj).NoPath;%length(s_info(jumpst,tj).instates);% time of jump is the same as tj
                    ind=ind+1;
                    % 						s_info(jumpst,tj).instates(ind)=s_i;
                    s_info(jumpst,tj).intime(ind)=t;
                    s_info(jumpst,tj).cost(ind)=s_info(s_i,t).cost(k)+sum((((codej*2-1).*code(t:tj-1))));%computing cost
                    s_info(jumpst,tj).path(ind).dec=[s_info(s_i,t).path(k).dec codej];
                    s_info(jumpst,tj).path(ind).ST=[s_info(s_i,t).path(k).ST jumpst];
                    s_info(jumpst,tj).NoPath=ind;
                    s_info(jumpst,tj).hasin=1;
                    chng(jumpst,tj)=1;
                else
                    ind=s_info(s_j,tj).NoPath;
                    ind=ind+1;
                    % 						s_info(s_j,tj).instates(ind)=s_i;
                    s_info(s_j,tj).intime(ind)=t;
                    s_info(s_j,tj).cost(ind)=s_info(s_i,t).cost(k)+sum((((codej*2-1).*code(t:tj-1))));%computing cost
                    s_info(s_j,tj).path(ind).dec=[s_info(s_i,t).path(k).dec codej];
                    s_info(s_j,tj).path(ind).ST=[s_info(s_i,t).path(k).ST s_j];
                    s_info(s_j,tj).NoPath=ind;
                    s_info(s_j,tj).hasin=1;
                    chng(s_j,tj)=1;
                end
            end
            %             s_info(s_j,tj)
        end
        %888888888888888888888888888888888888888888888
    end
end
% for state=1:sNo%show costs
%     s_info(state,t).cost
% end
highcost=[];k=1;
for state=1:sNo%show costs
    for in=1:s_info(state,t).NoPath
        if isequal(s_info(state,t).path(in).dec,R)
            disp(['state=',num2str(state),' time=',num2str(t),' input=',num2str(in)])
        end
        highcost(k).state=state;
        highcost(k).cost=s_info(state,t).cost(in);%sum([s_info(state,t).cost]/length(s_info(state,t).cost));
        %         highcost(k).instate=s_info(state,t).instates(in);
        highcost(k).intime=s_info(state,t).intime(in);
        highcost(k).path=s_info(state,t).path(in).dec;
        highcost(k).ST=s_info(state,t).path(in).ST;
        highcost(k).in=in;
        k=k+1;
    end
end

[sortedC,sortedind]=sort([highcost.cost],'descend');
highcost=highcost(sortedind);
highcost(L+1:end)=[];
