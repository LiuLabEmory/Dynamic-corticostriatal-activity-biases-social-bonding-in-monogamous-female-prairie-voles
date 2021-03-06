% script to generate cumulative summation plot for each of the animals

IDKeyFilePath = 'R:\LiuLab\People\Jim\Experiments\TwoChoiceOdorExposureExpt\AnimalIDandSide.xlsx';
% IDKeyFilePath = '/Volumes/ecas-research/LiuLab/People/Jim/Experiments/TwoChoiceOdorExposureExpt/AnimalIDandSide.xlsx';


[num,txt,~] = xlsread(IDKeyFilePath);

animID = cell2mat(txt(2:end,1));
animID_u = unique(animID,'rows');

BehStructNames = who('BehaviorStruct*');

for i = 1:length(animID_u)
    animID_temp = animID_u(i,:);
    animID_tf = cellfun(@(x) strcmp(x(16:17),animID_temp), BehStructNames);
    animID_i = find(animID_tf == 1);
    vidName = w{animID_i(1)};
    eval(['behStruct = ' vidName ]);
    
    if length(animID_i) > 1
        for vidNum = 2:length(animID_i) % concatenate all the videos
            eval(['behStruct_temp = ' w{animID_i(vidNum)}]);
            if max(behStruct_temp.Time_s) < 4000 % video is 30 min long
                behStruct_temp.Time_s = behStruct_temp.Time_s + (30*(vidNum-1)*60);
            else
                behStruct_temp.Time_s = behStruct_temp.Time_s + (52.5*(vidNum-1)*60);
            end
            behStruct.Behaviors = [behStruct.Behaviors; behStruct_temp.Behaviors];
            behStruct.Time_s = [behStruct.Time_s; behStruct_temp.Time_s];
            behStruct.Dur_s = [behStruct.Dur_s; behStruct_temp.Dur_s];
            % behStruct = [behStruct,behStruct_temp];
        end
    end
    
    sniff = [];
    side = [];
    

    sniff_i = behStruct.Behaviors == 1 | behStruct.Behaviors == 2;
    side_i = behStruct.Behaviors == 3 | behStruct.Behaviors == 4;
    
    sniff = behStruct.Dur_s(sniff_i);
    side = behStruct.Dur_s(side);
    
        % loop through all the behavior struct for that specific animal

            
end

        
