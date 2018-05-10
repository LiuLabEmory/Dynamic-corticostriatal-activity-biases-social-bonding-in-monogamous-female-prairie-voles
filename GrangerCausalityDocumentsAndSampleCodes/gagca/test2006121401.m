clear;
fnamestem='MEG10TailVz4RightV';
fnamein=[fnamestem '.mat'];
fnameout=[fnamestem 'afterAR.mat'];
load(fnamein);

[points,channelnum,n]=size(Data);
nobs0  = 300; 
Lmax0 = 15;
nspec = 300;
ncov = nobs0;
simux = 1;

for i=1:n
    for j=1:channelnum
        tmp=squeeze(Data(:,j,i));
        tmp=tmp-mean(tmp);
        Data(:,j,i)=tmp;
    end
end

for j=1:channelnum
    count=floor(n/5);
    for i=1:count
        tmp=squeeze(Data(:,j,(i-1)*5+1:i*5));
        tmp=reshape(tmp,points*5,1);
        tmp=lnrusa(tmp,600);
        tmp=reshape(tmp,points,5);
        Data(:,j,(i-1)*5+1:i*5)=tmp;
    end
end

for i=1:n
    for j=1:channelnum
        tmp=squeeze(Data(:,j,i));
        tmp=tmp-mean(tmp);
        tmp=tmp/std(tmp);
        Data(:,j,i)=tmp;
    end
end

for j=1:channelnum
    for i=1:points
        tmp=squeeze(Data(i,j,:));
        tmp=tmp-mean(tmp);
        tmp=tmp/std(tmp);
        Data(i,j,:)=tmp;
    end
end


x=squeeze(Data(601:900,3,:));
[points,n]=size(x);
for i=1:n
    y=x(:,i);
    Y = fft(y,300);
    if(i==1)
        Pyy = Y.* conj(Y) / 300;
    else
        Pyy(:,i) = Y.* conj(Y) /300;
    end
end
Pyym=mean(Pyy,2);

for taskstage=1:4
    OrderT=-1;
    for ARloop=1:n
        for i=1:channelnum
            v(i,1,:)=squeeze(Data((taskstage-1)*300+1:taskstage*300,i,ARloop));
        end

        %----------------------------------------------
        %Estimation
        %----------------------------------------------
        [pcsel, R0hat, logres, cic, pchat, psel, fsic] = ARselv(v,Lmax0);
        if(psel>OrderT)
            OrderT=psel;
        end
    end

    for ARloop=1:n
        for i=1:channelnum
            v(i,1,:)=squeeze(Data((taskstage-1)*300+1:taskstage*300,i,ARloop));
        end

        %----------------------------------------------
        %Estimation
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

    pcsel=mean(pcselAll,4);
    R0hat=mean(R0hatAll,3);
    %----------------------------------------------
    %spectra
    %----------------------------------------------
    [hsel gsel]= pc2specv(pcsel,R0hat,nspec);
    for l=1:channelnum
        for m=1:channelnum
            tmp=abs(squeeze(hsel(l,m,:)));
            tmp=tmp.*tmp;
            tmp=tmp./(squeeze(hsel(l,l,:)).*squeeze(hsel(m,m,:)));
            csel(l,m,:)=sqrt(tmp);
        end
    end
    
       
       
    if taskstage == 1
        Spectral=hsel;
        Granger=gsel;
        Order=OrderT;
        Coh=csel;
    else
        Spectral(:,:,:,taskstage)=hsel;
        Granger(:,:,:,taskstage)=gsel;
        Order(:,taskstage)=OrderT;
        Coh(:,:,:,taskstage)=csel;
    end
end

save(fnameout,'Pyym','Spectral','Granger','Order','Coh');

    

