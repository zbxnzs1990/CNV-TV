%2011/06/07
%jduan@tulane.edu
%using total variation penalized least squares to detect CNV from NGS data
%optimization problem is transformed to LASSO

clear all;
close all;
clc;


%%%%%%%%%%%%%%%%%%micro definitions%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_win=100;%bin size, depending on the inputed .hits file
TH=0.6;%threshold predifined that estimated from histogram
VAR=0.1;%noise variation predifined that estiamted from non CNV region
size_win2=1000;%bulk processing window size
D=toeplitz([1;zeros(size_win2-2,1)],[1,-1,zeros(1,size_win2-2)]);%first order derivative matrix
A=D'*inv(D*D');%pseudo inverse of D
B=toeplitz(ones(size_win2),[1,zeros(1,size_win2-1)]);


%%%%%%%%%%%%%%%%%%%%%load data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('ref.hits');
hits=fscanf(fid,'%d \n',[1,inf]);
% fid=fopen('J:\share_linux_suse\hg18chr21.hits');
% hits=fscanf(fid,'%c %c %c %d %d \n',[5,inf]);
fclose(fid);
loc=hits(1,:);


%%%%%%%%%%%%converting reference hits to read depth if necessary%%%%%%%%%%%%%%
tmp=zeros(max(loc),1);
for i=1:length(loc)
    tmp(loc(i))=tmp(loc(i))+1;
end

rd_ref=zeros(floor(max(loc)/size_win),1);
for i=1:floor(max(loc)/size_win)
    rd_ref(i)=sum(tmp((i-1)*size_win+1:i*size_win));
end


%%%%%%%%%%%%%%%converting test hits to read depth%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('test.hits');
hits=fscanf(fid,'%d \n',[1,inf]);
%fid=fopen('J:\share_linux_suse\NA19240_HiSeq_100_chr21.hits');
%hits=fscanf(fid,'%c %c %c %d %d \n',[5,inf]);
fclose(fid);
loc=hits(1,:);

tmp=zeros(max(loc),1);
for i=1:length(loc)
    tmp(loc(i))=tmp(loc(i))+1;
end

rd_test=zeros(floor(max(loc)/size_win),1);
for i=1:floor(max(loc)/size_win)
    rd_test(i)=sum(tmp((i-1)*size_win+1:i*size_win));
end


%%%%%%%%%%%%%%%%%%%%%%%filter out no data regions%%%%%%%%%%%%%%%%%%%%%%
len=min([length(rd_ref),length(rd_test)]);
rd_ratio=rd_test(1:len)./rd_ref(1:len);

idx=find((~isnan(rd_ratio)&(~isinf(rd_ratio))));
rd_ratio=rd_ratio/mean(rd_ratio(idx));
%idx=find(isinf(rd_ratio));
% rd_ratio(idx)=1;
% idx=find(isnan(rd_ratio));
% rd_ratio(idx)=1;
% idx=find(rd_ratio==0);
% rd_ratio(idx)=1;
%rd=log2(rd_ratio);
rd=rd_ratio;
rd=rd/mean(rd);


%%%%%%%%%%%%%%%%%%%piece wise smoothing%%%%%%%%%%%%%%%%%%%%%%
idx_bp=[];
for i=1:floor(length(rd)/size_win2)
%    i
    rd1=rd((i-1)*size_win2+1:i*size_win2);
    u=SolveLasso_Flops(A,A*D*rd1,'lasso',100,1);%call standard LASSO solver
    figx=[];
    figy=[];
    for j=1:size(u,2)%re least square
        idx=find(u(:,j));
        C=B(:,[1;idx]);
        figy(j)=rd1'*rd1-rd1'*C*pinv(C'*C)*C'*rd1;
        figy(j)=figy(j)/VAR+log(size_win2)*length(idx);%calculate SIC
        figx(j)=length(idx);
    end
    
    [NA,idx_opt]=min(figy);
    bps=find(u(:,idx_opt))';%store breakpoint locations
    if ~isempty(bps)
        idx_bp=[idx_bp,(i-1)*size_win2+bps];
    end
    
    idx=find(u(:,idx_opt));
    C=B(:,[1;idx]);
    yy=C*pinv(C)*rd1;

    
    %%%%%%%%%%%%%%%%%draw SIC curve for model selection%%%%%%%%%%%%%%%%
    figure(1);
    plot(figx,figy,'.');
    hold on;
    plot(figx(idx_opt),figy(idx_opt),'ro','linewidth',1.5);
    xlabel('k');
    ylabel('SIC');
    drawnow;
    hold off;
 
    
    %%%%%%%%%%%%%%%%%draw smoothed signal%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2);
    plot((i*size_win2+[-(size_win2-1):0])*size_win,rd1,'k.','markersize',4);
    hold on;
    plot((i*size_win2+[-(size_win2-1):0])*size_win,yy,'r');
    plot((i*size_win2+[-(size_win2-1):0])*size_win,rd1+A*(u(:,idx_opt)-D*rd1),'b');
    xlabel('chromosome position');
    ylabel('read depth');
    drawnow;
    hold off;
    pause;
end

clear A B D;


%%%%%%%%%%%%%%%%%%insert boders as potetion break points%%%%%%%%%%%%%%%%
idx_bp=[idx_bp,size_win2:size_win2:length(rd)];
idx_bp=sort(idx_bp);
idx_bp=[0,idx_bp,length(rd)];
idx_bp=unique(idx_bp);


%%%%%%%%%%%%%%%%%filter out segments with small size%%%%%%%%%%%%%%%%%% 
for i=length(idx_bp):-1:2
    if (abs(idx_bp(i)-idx_bp(i-1))<10)
        idx_bp(i)=[];
    end
end


%%%%%%%%%%%%%%%%%piece wise smoothing%%%%%%%%%%%%%%%%%%%%%%%%%%
dd=[];
for i=1:length(idx_bp)-1
    dd(idx_bp(i)+1:idx_bp(i+1))=mean(rd(idx_bp(i)+1:idx_bp(i+1)));
end


%%%%%%%%%%%%%%%%%%%%%summarize calls%%%%%%%%%%%%%%%%%%%%%%%%%%
amp=unique(dd);
idx_amp=find(abs(amp-1)>TH);
if ~isempty(idx_amp)
    for i=1:length(idx_amp)
        loc_cnv=find(dd==amp(idx_amp(i)));
        pos_begin=(loc_cnv(1)-1)*size_win+1;
        pos_end=loc_cnv(end)*size_win;
%        cpn=amp(idx_amp(i))/mean(rd);
        cpn=amp(idx_amp(i));
        out(i,:)=[pos_begin,pos_end,cpn];
    end
    
    [NA,idx]=sort(out(:,1));
    out1=out(idx,:);
    out=[];
    
    lab(1)=1;
    for i=1:size(out1,1)-1
        if (out1(i,2)+1)==out1(i+1,1)
            lab(i+1)=lab(i);
        else
            lab(i+1)=lab(i)+1;
        end
    end
    
    for i=1:lab(end)
        idx=find(lab==i);
        out(i,1)=out1(idx(1),1);
        out(i,2)=out1(idx(end),2);
        out(i,3)=sum((out1(idx,2)-out1(idx,1)+1).*out1(idx,3))/(out(i,2)-out(i,1));
    end
    
    
%%%%%%%%%%%%%%%%%%plot calls%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
plot(rd);
hold on;
xlabel('position (100bp)','fontsize',12);
ylabel('read depth','fontsize',12);
plot(dd,'r');
plot(1:length(dd),1+TH,'g');
plot(1:length(dd),1-TH,'g');

for i=1:size(out,1)
    idx=((out(i,1)-1)/size_win+1):(out(i,2)/size_win);
    plot(idx,dd(idx),'r','linewidth',5);
end


%%%%%%%%%%%%%%%%%%%%%print calls into file%%%%%%%%%%%%%%%%%%%    
    fid=fopen('J:\share_linux_suse\out1.txt','w');
    fprintf(fid,'start\tend\tcopy number\n');
    for i=1:size(out,1)
        fprintf(fid,'%d\t%d\t%f\n',out(i,:));
    end
    fclose(fid);
end
