% GrangerZxliuARMovinggg2

%% Input parameters

LFPdirect = 'C:\Georgia\PlexonData\lfpsnewFilt\';
outdir = 'C:\Georgia\PlexonData\Granger\cueon\contraout\timecourse\win150step10\neg200_1000\';

maxsteps=121;
windowstep=10;

%the column number of each task related information in the data matrix
colorchange=0;
trialnumcol=1;
starttime=2;
endtime=3;
stimulusontime=4;
targetchange=5;
correctornot=6;
cueontime=7;
targcolor=8;
targx=9;
targy=10;
stimtocue=11;
stimtotargchange=12; % time from stimulusontime to stimuluschangetime
distractor1change=13;
distractor2change=14;
distractor1x=15;
distractor1y=16;
distractor1color=17;
distractor2x=18;
distractor2y=19;
distractor2color=20;
dimcondition=21;
cuetotargchange=22;
%%% the following information is for memsacc task
memsaccstimon=23;
memsaccgocue=24;
memsaccx=25;
memsaccy=26;

%% Load files

cd('C:\Georgia\PlexonData\mua_coh\lfp_lfp\cueon\');


lfplfpinfo = listfromfile('lfplfpcohlistFEFV4cue_popaver.txt');

lfplfpnumber=length(lfplfpinfo);


fid = fopen('lfplfpcohlistFEFV4cue_popaver.txt','r');
if(fid == -1)
    disp('cannot open file');
    return
end

for pairnumber=1:lfplfpnumber
    linenumber=pairnumber;
    currline=lfplfpinfo{pairnumber};
    disp(currline);
    infopositions=findstr(currline,':');
    LFP1fname=currline(infopositions(1)+1:infopositions(2)-1);
    LFP2fname=currline(infopositions(2)+1:infopositions(3)-1);
    RFX1=str2num(currline(infopositions(3)+1:infopositions(4)-1));
    RFY1=str2num(currline(infopositions(4)+1:infopositions(5)-1));

    outx1=str2num(currline(infopositions(5)+1:infopositions(6)-1));
    outy1=str2num(currline(infopositions(6)+1:infopositions(7)-1));
    outx2=str2num(currline(infopositions(9)+1:infopositions(10)-1));
    outy2=str2num(currline(infopositions(10)+1:infopositions(11)-1));



    %%% load LFP1 data
    LFP1filename=strcat(LFPdirect,LFP1fname);
    LFP2filename=strcat(LFPdirect,LFP2fname);


    cd(LFPdirect);
    if exist(LFP1filename)&& exist(LFP2filename)
        load(LFP1filename); % get LFP data
        load(dataforLFP);     % get trial information (on trial events) associated with LFPs
        data=informationdata1; 
        clear informationdata1;
        clear dataforLFP;
        clear AD_15;
        clear AD_16;
        clear AD15_ts;
        clear AD16_ts;
        clear adfreq;

        LFP_startind=round(initialtime*1000);

        clear fn;
        clear initialtime;


        totallength=length(LFPsignal);
        count=ceil(totallength/(1000.0*10));
        for j=1:count
            startind=(j-1)*1000*10+1;
            endind=j*1000*10;
            if(endind>totallength)
                endind=totallength;
            end
            tmp=LFPsignal(startind:endind);
            tmp=lnrusa4(tmp,1000);
            LFPsignal(startind:endind)=tmp;
        end


        LFPsignalA=LFPsignal';

        clear LFPsignal 



        %%%% load LFP2 data

        % disp(LFP2filename);
        load(LFP2filename); % get LFP data
        load(dataforLFP);     % get trial information (on trial events) associated with LFPs
        data=informationdata1; % 
        clear informationdata1;
        clear dataforLFP;
        clear AD_15;
        clear AD_16;
        clear AD15_ts;
        clear AD16_ts;
        clear adfreq;

        LFP_startind=round(initialtime*1000);
        clear fn;
        clear initialtime;


        totallength=length(LFPsignal);
        count=ceil(totallength/(1000.0*10));
        for j=1:count
            startind=(j-1)*1000*10+1;
            endind=j*1000*10;
            if(endind>totallength)
                endind=totallength;
            end
            tmp=LFPsignal(startind:endind);
            tmp=lnrusa(tmp,1000);
            LFPsignal(startind:endind)=tmp;
        end


        LFPsignalB=LFPsignal';

        clear LFPsignal 


        %% select attendin trials and find start and end indices

            AttendInTypeATrials=find(data(:,correctornot)==200 & data(:,targx)==RFX1 & data(:,targy)==RFY1 & data(:,stimtocue) <0);   % attendin trials ,cue on before grating stimulus onset
            TrialNum=length(AttendInTypeATrials);
            kk=0;
            for j=1:TrialNum
                cueindex=data(AttendInTypeATrials(j),cueontime);
                cueindex=round(cueindex*1000)-LFP_startind+1;
                endindex=data(AttendInTypeATrials(j),targetchange); %coherence calculation end just before stimulus changed
                endindex=round(endindex*1000)-LFP_startind+1;
                if endindex-stiindex>=1000
                    kk=kk+1;
                    TypeAStartin(kk,1)=stiindex-275;
                end
            end


        %% select attendout trials and find start and end indices

            AttendOutTypeATrials=find(data(:,correctornot)==200 & (data(:,targx)==outx1 & data(:,targy)==outy1) & data(:,stimtocue) < 0);   % attendout trials ,cue on before grating stimulus onset


            TrialNum=length(AttendOutTypeATrials);
            
            kk=0;
            for j=1:TrialNum
                cueindex=data(AttendOutTypeATrials(j),cueontime);
                cueindex=round(cueindex*1000)-LFP_startind+1;

                endindex=data(AttendOutTypeATrials(j),targetchange); %coherence calculation end just before stimulus changed
                endindex=round(endindex*1000)-LFP_startind+1;
                
                if endindex-stiindex>=1000
                kk=kk+1;
                TypeAStartout(kk,1)=stiindex-275;
                end
            end

            Startout=TypeAStartout;
            Startin=TypeAStartin;


        
%         %% select equal number of trials for attendin and attendout
%         
% 
% 
%         if length(Startout)> length(Startin)
%             index=randperm(length(Startout));
%             index=index';
%             selecttrial=find(index(:,1)<=length(Startin));
%             Startout=Startout(selecttrial,1);
%             clear index;
%             clear selecttrial;
%         elseif length(Startout)<= length(Startin)
%             TrialNum=length(Startout);
%             Startin=Startin(1:TrialNum,:);
%             clear index;
%             clear selecttrial;
%         end
%         

            TypeAStart=[Startin Startout];
  

        %% Type A calculation

            [TrialNum,att]=size(TypeAStart);

            stage=maxsteps; % total window steps
            for ia=1:stage %do entire analysis separately for each window and attention condition
                for ja=1:att
                    
                    
                    
                    Data=zeros(150,2,TrialNum);
                    for ka=1:TrialNum
                        Data(:,1,ka)=LFPsignalA(TypeAStart(ka,ja)+(ia-1)*windowstep:TypeAStart(ka,ja)+(ia-1)*windowstep+149);%150 datapts, LFP signal is trials x condition 
                        Data(:,2,ka)=LFPsignalB(TypeAStart(ka,ja)+(ia-1)*windowstep:TypeAStart(ka,ja)+(ia-1)*windowstep+149);
                    end

                    %----------------------------------------------
                    %Preprocessing
                    %----------------------------------------------

                    [points,channelnum,n]=size(Data);

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
                    for ARloop=1:n
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

                    for ARloop=1:n
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

                    pcsel=mean(pcselAll,4);
                    R0hat=mean(R0hatAll,3);

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

                    TypeACoh(ia,ja,:)=csel;
                    TypeAGrangerDown(ia,ja,:)=gtopdownsel;
                    TypeAGrangerUp(ia,ja,:)=gbottomupsel;
                    TypeAOrder(ia,ja,:)=OrderT;
                end
            end


        clear OrderT gbottomupsel gtopdownsel gseltmp hsel pcsel R0hat nspec pcsellAll R0hatAll v logres cic pchat psel fsic ARloop Data

        %%


        clear data;
        clear LFPsignalA;
        clear LFPsignalB;
        clear TypeAStartin TypeAStartout TypeAStart Startin Startout;
        clear TypeAEnd;

        fname=[outdir LFP1fname(1,1:end-4)  '_' LFP2fname(1,end-9:end)];

            save(fname, 'TypeACoh', 'TypeAGrangerDown', 'TypeAGrangerUp', 'TypeAOrder')

        clear TypeACoh TypeAGrangerDown TypeAGrangerUp TypeAOrder tmp LFP1signal LFP2signal AD15 AD16...
            LFP1filename LFP2filename LFP1fname LFP2fname

    end
end

