% N Killian, edited down, circa 11-2010


%----------------------------------------------
%Preprocessing
%----------------------------------------------

[points,channelnum,n]=size(Data);%time samples x number of channels x trials *this is correct*

Lmax0 = 30;  % max order to try in order selection
nspec = 500;

% subtract mean of each segment in a trial and divide by std
for i=1:n
    for j=1:channelnum
        tmp=squeeze(Data(:,j,i));
        tmp=tmp-mean(tmp);
        tmp=tmp/std(tmp);
        Data(:,j,i)=tmp;
    end
end

% subtract mean of all trials per ms and divide by std
for j=1:channelnum
    for i=1:points
        tmp=squeeze(Data(i,j,:));
        tmp=tmp-mean(tmp);
        tmp=tmp/std(tmp);
        Data(i,j,:)=tmp;
    end
end

OrderT=-1;
for ARloop=1:n %n = number of trials
    for i=1:channelnum
        v(i,1,:)=squeeze(Data(:,i,ARloop));% v is chans x 1 x samples (1 trial)
    end
    
    %----------------------------------------------
    %Estimation/Selection of model order
    %----------------------------------------------
    [pcsel, R0hat, logres, cic, pchat, psel, fsic] = ARselv(v,Lmax0);
    if(psel>OrderT)
        OrderT=psel;
    end
end

for ARloop=1:n% all trials
    for i=1:channelnum
        v(i,1,:)=squeeze(Data(:,i,ARloop));
    end
    
    %----------------------------------------------
    %Estimation of model's parameters for the selected order
    %----------------------------------------------
    [pcsel, R0hat, logres, cic, pchat, psel, fsic] = ARselvFixedOrder(v,OrderT);
    
    
    if ARloop==1
        pcselAll=pcsel;
        R0hatAll=R0hat;
    else
        pcselAll(:,:,:,ARloop)=pcsel;
        R0hatAll(:,:,ARloop)=R0hat;
    end
end

pcsel=mean(pcselAll,4);%avg across trials, channel x channel x time?
R0hat=mean(R0hatAll,3);% avg across trials, trial x channel

%----------------------------------------------
%spectra
%----------------------------------------------
[hsel gsel]= pc2specv(pcsel,R0hat,nspec);
tmp=abs(squeeze(hsel(1,2,:)));
tmp=tmp.*tmp;
tmp=tmp./(squeeze(hsel(1,1,:)).*squeeze(hsel(2,2,:)));
csel=sqrt(tmp);
gtopdownsel=squeeze(gsel(1,2,:));
gbottomupsel=squeeze(gsel(2,1,:));

TypeACoh(ia,ja,:)=csel;%step x condition x freqs
TypeAGrangerDown(ia,ja,:)=gtopdownsel;
TypeAGrangerUp(ia,ja,:)=gbottomupsel;
TypeAOrder(ia,ja,:)=OrderT;


