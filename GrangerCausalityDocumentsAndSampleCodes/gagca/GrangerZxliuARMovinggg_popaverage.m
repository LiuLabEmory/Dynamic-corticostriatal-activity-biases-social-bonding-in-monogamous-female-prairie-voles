% GrangerZxliuARMovinggg_popaverage


zeropoints=input('which zero time: [4: stimulusontime] ; [5: targetchange] ; [7: cueon];');
bin=input('use window average for statistics: yes: 1 ; no: 0 ;');
stat=input('plot only datapoints with statistically significant attn effect: 1 ; plot all coherence differences: 0 ; ');
Grangerdirect = 'C:\Georgia\PlexonData\Granger\cueon\contraout\timecourse\win150step10\neg200_1000\';
cd('C:\Georgia\PlexonData\mua_coh\lfp_lfp\cueon\');

maxfreq=150;
maxsteps=121;
p=0.05;
windowsize=1;
consbin=10;
slidingwindowsize=150;
startpoint=275;
step=10;

lfplfpinfo = listfromfile('lfplfpcohlistFEFV4cue_popaver.txt');
lfplfpnumber=length(lfplfpinfo);

fid = fopen('lfplfpcohlistFEFV4cue_popaver.txt','r');
if(fid == -1)
    disp('cannot open file');
    return
end

k=0;
for Grangerfile=1:lfplfpnumber
    currline=lfplfpinfo{Grangerfile};
    disp(currline);
    infopositions=findstr(currline,':');
    Grangerfname=strcat(currline(infopositions(1)+1:infopositions(2)-5),'_',currline(infopositions(3)-10:infopositions(3)-1));

    Grangerfilename=strcat(Grangerdirect,Grangerfname);

    if exist(Grangerfilename)
        k=k+1;
        load(Grangerfilename);

            Coh=TypeACoh;
            GrangerDown=TypeAGrangerDown;
            GrangerUp=TypeAGrangerUp;
            Order=TypeAOrder;

        if Grangerfile==1
            CohInAll=squeeze(Coh(1:maxsteps,1,1:maxfreq));
            GrangerDownInAll=squeeze(GrangerDown(1:maxsteps,1,1:maxfreq));
            GrangerUpInAll=squeeze(GrangerUp(1:maxsteps,1,1:maxfreq));
            OrderInAll=Order(1:maxsteps,1);
            CohOutAll=squeeze(Coh(1:maxsteps,2,1:maxfreq));
            GrangerDownOutAll=squeeze(GrangerDown(1:maxsteps,2,1:maxfreq));
            GrangerUpOutAll=squeeze(GrangerUp(1:maxsteps,2,1:maxfreq));
            OrderOutAll=Order(1:maxsteps,2);

        elseif Grangerfile>1
            CohInAll(:,:,k)=squeeze(Coh(1:maxsteps,1,1:maxfreq));
            GrangerDownInAll(:,:,k)=squeeze(GrangerDown(1:maxsteps,1,1:maxfreq));
            GrangerUpInAll(:,:,k)=squeeze(GrangerUp(1:maxsteps,1,1:maxfreq));
            OrderInAll(:,k)=Order(1:maxsteps,1);
            CohOutAll(:,:,k)=squeeze(Coh(1:maxsteps,2,1:maxfreq));
            GrangerDownOutAll(:,:,k)=squeeze(GrangerDown(1:maxsteps,2,1:maxfreq));
            GrangerUpOutAll(:,:,k)=squeeze(GrangerUp(1:maxsteps,2,1:maxfreq));
            OrderOutAll(:,k)=Order(1:maxsteps,2);

            clear Coh GrangerDown GrangerUp Order

        end

    end


end

[steps  freq n]=size(CohInAll);
CohInAll=abs(CohInAll);
CohInMean=mean(CohInAll,3);
CohOutAll=abs(CohOutAll);
CohOutMean=mean(CohOutAll,3);

CohAtt=CohInMean-CohOutMean;


if stat==1
    for i=1:steps
        for j=1:freq
            x=squeeze(CohInAll(i,j,:));
            y=squeeze(CohOutAll(i,j,:));
            CohP(i,j)=ttest(x,y,p);
        end
    end
    CohAtt=CohAtt.*CohP;
end

freq=1:maxfreq;
subplot(3,1,1)
tmp=CohAtt(:,freq)';
[freqs steps]=size(tmp);
time=1:steps;
time=time*step-startpoint+slidingwindowsize/2;
imagesc(time,freq,tmp);
axis xy; colorbar;

[steps  freq n]=size(GrangerDownInAll);
GrangerDownInAll=abs(GrangerDownInAll);
GrangerDownInMean=mean(GrangerDownInAll,3);
GrangerDownOutAll=abs(GrangerDownOutAll);
GrangerDownOutMean=mean(GrangerDownOutAll,3);

GrangerDownAtt=GrangerDownInMean-GrangerDownOutMean;


if stat==1
    for i=1:steps
        for j=1:freq
            x=squeeze(GrangerDownInAll(i,j,:));
            y=squeeze(GrangerDownOutAll(i,j,:));
            GrangerDownP(i,j)=ttest(x,y,p);
        end
    end
    GrangerDownAtt=GrangerDownAtt.*GrangerDownP;
end

freq=1:maxfreq;
subplot(3,1,2)
tmp=GrangerDownAtt(:,freq)';
[freqs steps]=size(tmp);
time=1:steps;
time=time*step-startpoint+slidingwindowsize/2;
imagesc(time,freq,tmp);
axis xy; colorbar;

[steps freq n]=size(GrangerUpInAll);
GrangerUpInAll=abs(GrangerUpInAll);
GrangerUpInMean=mean(GrangerUpInAll,3);
GrangerUpOutAll=abs(GrangerUpOutAll);
GrangerUpOutMean=mean(GrangerUpOutAll,3);

GrangerUpAtt=GrangerUpInMean-GrangerUpOutMean;


if stat==1
    for i=1:steps
        for j=1:freq
            x=squeeze(GrangerUpInAll(i,j,:));
            y=squeeze(GrangerUpOutAll(i,j,:));
            GrangerUpP(i,j)=ttest(x,y,p);
        end
    end
    GrangerUpAtt=GrangerUpAtt.*GrangerUpP;
end

freq=1:maxfreq;
subplot(3,1,3)
tmp=GrangerUpAtt(:,freq)';
[freqs steps]=size(tmp);
time=1:steps;
time=time*step-startpoint+slidingwindowsize/2;
imagesc(time,freq,tmp);
axis xy; colorbar;

figure
subplot(2,2,1)
imagesc(time,freq,GrangerDownInMean');
axis xy;colorbar;
subplot(2,2,2)
imagesc(time,freq,GrangerDownOutMean');
axis xy;colorbar;

subplot(2,2,3)
imagesc(time,freq,GrangerUpInMean');
axis xy;colorbar;
subplot(2,2,4)
imagesc(time,freq,GrangerUpOutMean');
axis xy;colorbar;


%% Normalize to gamma baseline

[steps  freq n]=size(GrangerUpInAll);

for i=1:n
    GrangerUpInAllnormal(:,:,i)= GrangerUpInAll(:,:,i)/mean([mean(mean(GrangerUpInAll(1:20,40:60,i),2)) mean(mean(GrangerUpOutAll(1:20,40:60,i),2))]);
    GrangerUpOutAllnormal(:,:,i)= GrangerUpOutAll(:,:,i)/mean([mean(mean(GrangerUpInAll(1:20,40:60,i),2)) mean(mean(GrangerUpOutAll(1:20,40:60,i),2))]);
    GrangerDownInAllnormal(:,:,i)= GrangerDownInAll(:,:,i)/mean([mean(mean(GrangerDownInAll(1:20,40:60,i),2)) mean(mean(GrangerDownOutAll(1:20,40:60,i),2))]);
    GrangerDownOutAllnormal(:,:,i)= GrangerDownOutAll(:,:,i)/mean([mean(mean(GrangerDownInAll(1:20,40:60,i),2)) mean(mean(GrangerDownOutAll(1:20,40:60,i),2))]);
    for ii=1:steps
        GrangerUpAttdiffnorm(ii,i)=mean(GrangerUpInAllnormal(ii,40:60,i))-mean(GrangerUpOutAllnormal(ii,40:60,i));
        GrangerDownAttdiffnorm(ii,i)=mean(GrangerDownInAllnormal(ii,40:60,i))-mean(GrangerDownOutAllnormal(ii,40:60,i));
        [h_p(ii),p_p(ii)]=ttest(GrangerUpAttdiffnorm(ii,:),GrangerDownAttdiffnorm(ii,:));
    end

end

GrangerUpInMeannormal=mean(GrangerUpInAllnormal,3);
GrangerUpOutMeannormal=mean(GrangerUpOutAllnormal,3);
GrangerDownInMeannormal=mean(GrangerDownInAllnormal,3);
GrangerDownOutMeannormal=mean(GrangerDownOutAllnormal,3);

GrangerUpAttnormal=GrangerUpInMeannormal-GrangerUpOutMeannormal;
GrangerDownAttnormal=GrangerDownInMeannormal-GrangerDownOutMeannormal;

if stat==1
    for i=1:steps
        for j=1:freq
            x=squeeze(GrangerUpInAllnormal(i,j,:));
            y=squeeze(GrangerUpOutAllnormal(i,j,:));
            GrangerUpPnormal(i,j)=ttest(x,y,p);
        end
    end
    GrangerUpAttnormal=GrangerUpAttnormal.*GrangerUpPnormal;
end

figure
freq=1:maxfreq;
subplot(2,1,2)
tmp=GrangerUpAttnormal(:,freq)';
[freqs steps]=size(tmp);
time=1:steps;
time=time*step-startpoint+slidingwindowsize/2;
imagesc(time,freq,tmp);
axis xy; colorbar;

[steps  freq n]=size(GrangerUpInAll);

if stat==1
    for i=1:steps
        for j=1:freq
            x=squeeze(GrangerDownInAllnormal(i,j,:));
            y=squeeze(GrangerDownOutAllnormal(i,j,:));
            GrangerDownPnormal(i,j)=ttest(x,y,p);
        end
    end
    GrangerDownAttnormal=GrangerDownAttnormal.*GrangerDownPnormal;
end

freq=1:maxfreq;
subplot(2,1,1)
tmp=GrangerDownAttnormal(:,freq)';
[freqs steps]=size(tmp);
time=1:steps;
time=time*step-startpoint+slidingwindowsize/2;
imagesc(time,freq,tmp);
axis xy; colorbar;

title('Normalized to gamma baseline');


%%

%----------------------------------

GrangerUpOutgamma=mean(GrangerUpOutAll(:,40:60,:),2);
GrangerUpIngamma=mean(GrangerUpInAll(:,40:60,:),2);
sterrgrangerupoutgamma=std(squeeze(GrangerUpOutgamma(:,:)),0,2)/sqrt(k);
sterrgrangerupingamma=std(squeeze(GrangerUpIngamma(:,:)),0,2)/sqrt(k);

for i=1:size(GrangerUpOutgamma,3)
    GrangerUpOutgammanormal(:,i)=squeeze(GrangerUpOutgamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(1:20,1,i)),1) mean(squeeze(GrangerUpIngamma(1:20,1,i)),1)]);
    GrangerUpIngammanormal(:,i)=squeeze(GrangerUpIngamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(1:20,1,i)),1) mean(squeeze(GrangerUpIngamma(1:20,1,i)),1)]);
end
sterrgrangerupoutgammanormal=std(GrangerUpOutgammanormal,0,2)/sqrt(k);
sterrgrangerupingammanormal=std(GrangerUpIngammanormal,0,2)/sqrt(k);

    GrangerUp=(mean(squeeze(GrangerUpIngamma),2)-mean(squeeze(GrangerUpOutgamma),2));

    GrangerUpnormal=(mean(GrangerUpIngammanormal,2)-mean(GrangerUpOutgammanormal,2));


figure
subplot(2,1,2)
plot(time,mean(squeeze(GrangerUpIngamma),2),'r')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpIngamma),2)+sterrgrangerupingamma)' fliplr((mean(squeeze(GrangerUpIngamma),2)-sterrgrangerupingamma)')],'r', 'FaceAlpha',.3, 'EdgeColor','r','EdgeAlpha',.1);
plot(time, mean(squeeze(GrangerUpOutgamma),2),'b')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpOutgamma),2)+sterrgrangerupoutgamma)' fliplr((mean(squeeze(GrangerUpOutgamma),2)-sterrgrangerupoutgamma)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);


GrangerDownOutgamma=mean(GrangerDownOutAll(:,40:60,:),2);
GrangerDownIngamma=mean(GrangerDownInAll(:,40:60,:),2);
sterrgrangerdowningamma=std(squeeze(GrangerDownIngamma(:,:)),0,2)/sqrt(k);
sterrgrangerdownoutgamma=std(squeeze(GrangerDownOutgamma(:,:)),0,2)/sqrt(k);

for i=1:size(GrangerDownOutgamma,3)
    GrangerDownOutgammanormal(:,i)=squeeze(GrangerDownOutgamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(1:20,1,i)),1) mean(squeeze(GrangerDownIngamma(1:20,1,i)),1)]);
    GrangerDownIngammanormal(:,i)=squeeze(GrangerDownIngamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(1:20,1,i)),1) mean(squeeze(GrangerDownIngamma(1:20,1,i)),1)]);
end
sterrgrangerdownoutgammanormal=std(GrangerDownOutgammanormal,0,2)/sqrt(k);
sterrgrangerdowningammanormal=std(GrangerDownIngammanormal,0,2)/sqrt(k);


    GrangerDown=(mean(squeeze(GrangerDownIngamma),2)-mean(squeeze(GrangerDownOutgamma),2));

    GrangerDownnormal=(mean(GrangerDownIngammanormal,2)-mean(GrangerDownOutgammanormal,2));

subplot(2,1,1)
plot(time,mean(squeeze(GrangerDownIngamma),2),'r')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerDownIngamma),2)+sterrgrangerdowningamma)' fliplr((mean(squeeze(GrangerDownIngamma),2)-sterrgrangerdowningamma)')],'r', 'FaceAlpha',.3, 'EdgeColor','r','EdgeAlpha',.1);
hold on
plot(time, mean(squeeze(GrangerDownOutgamma),2),'b')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerDownOutgamma),2)+sterrgrangerdownoutgamma)' fliplr((mean(squeeze(GrangerDownOutgamma),2)-sterrgrangerdownoutgamma)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);


figure
subplot(2,1,2)
plot(time,mean(GrangerUpIngammanormal,2),'r')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpIngammanormal),2)+sterrgrangerupingammanormal)' fliplr((mean(squeeze(GrangerUpIngammanormal),2)-sterrgrangerupingammanormal)')],'r', 'FaceAlpha',.3, 'EdgeColor','r','EdgeAlpha',.1);
hold on
plot(time, mean(GrangerUpOutgammanormal,2),'b')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpOutgammanormal),2)+sterrgrangerupoutgammanormal)' fliplr((mean(squeeze(GrangerUpOutgammanormal),2)-sterrgrangerupoutgammanormal)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);


subplot(2,1,1)
plot(time,mean(GrangerDownIngammanormal,2),'r')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerDownIngammanormal),2)+sterrgrangerdowningammanormal)' fliplr((mean(squeeze(GrangerDownIngammanormal),2)-sterrgrangerdowningammanormal)')],'r', 'FaceAlpha',.3, 'EdgeColor','r','EdgeAlpha',.1);
hold on
plot(time, mean(GrangerDownOutgammanormal,2),'b')
hold on
fill([time fliplr(time)] ,[(mean(squeeze(GrangerDownOutgammanormal),2)+sterrgrangerdownoutgammanormal)' fliplr((mean(squeeze(GrangerDownOutgammanormal),2)-sterrgrangerdownoutgammanormal)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);


title('Normalized to baseline');


