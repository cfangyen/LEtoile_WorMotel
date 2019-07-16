% Matt Churgin and Chris Fang-Yen
% Last edited 11/29/17
%% Load Data
clear all
close all

[srcfilename pname]=uigetfile('*.mat','MultiSelect', 'on','Select the .mat files');
cd(pname);
load(srcfilename);

group1 = [1 2 3 4 9 10 11 12 17 18 19 20 25 26 27 28 33 34 35 36 41 42 43 44 ];  % TRAINING

group2 = [5 6 7 8 13 14 15 16 21 22 23 24 29 30 31 32 37 38 39 40 45 46 47 48 ];   % CONTROL


time_between_frames = 3; % time between images in seconds

ActValcrop1 = ActVal(group1, :);
ActValcrop2 = ActVal(group2, :);

%%

numworms1 = length(group1);
numworms2 = length(group2);

quiescence_threshold = 1;

numtimepoints = size(ActVal,2);
fontsize2 = 15;

figure(1); 
subplot(211);
imagesc(ActValcrop1);
title('Group 1 activity')
colorbar;
xlabel('Frame #', 'FontSize', fontsize2)
ylabel('Worm', 'FontSize', fontsize2)
subplot(212);
imagesc(ActValcrop2);
title('Group 2 activity')
colorbar;
xlabel('Frame #', 'FontSize', fontsize2)
ylabel('Worm', 'FontSize', fontsize2)

ActValcrop1Q = (ActValcrop1 >= quiescence_threshold);
ActValcrop2Q = (ActValcrop2 >=  quiescence_threshold);
figure(2);
subplot(211);
imagesc(ActValcrop1Q);
subplot(212);
imagesc(ActValcrop2Q);
ActValcrop1totalQ = sum(1-ActValcrop1Q,2);
ActValcrop2totalQ = sum(1-ActValcrop2Q,2);

[ ActValcrop1Qsorted, index1] = sortrows([ActValcrop1totalQ ActValcrop1Q ]);
[ ActValcrop2Qsorted, index2] = sortrows([ActValcrop2totalQ ActValcrop2Q ]);

ActValcrop1totalQsorted = ActValcrop1Qsorted(:,1)* time_between_frames/60;
ActValcrop2totalQsorted = ActValcrop2Qsorted(:,1)* time_between_frames/60;
ActValcrop1Qsorted = ActValcrop1Qsorted(:,2:end);
ActValcrop2Qsorted = ActValcrop2Qsorted(:,2:end);

fontsize1 = 25;
fontsize2  = 12;

figure(3);clf;
subplot(221);
imagesc(ActValcrop1Qsorted);
xlabel('Time (min)', 'FontSize', fontsize1)
ax = gca;
ax.XTick = 0:400:numtimepoints;
ax.XTickLabel = [0:400:numtimepoints] * time_between_frames/60;
ylabel('Worm Number', 'FontSize', fontsize1)
ax = gca;
ax.YTick = 1:numworms1;
% ax.YTickLabel = ActValcrop1totalQsorted* time_between_frames/60;
ax.YTickLabel = group1(index1); % Comment this out for sequential labeling
%title(PathName, 'Interpreter', 'None', 'FontSize', fontsize1);
title('Group 1 (Quiescent = dark)');
set(gca,'FontSize',fontsize2);
% comment next three lines out to remove quiescence labels on right side of heat map
for j=1:numworms1
    text(numtimepoints*1.01,j, num2str(ActValcrop1totalQsorted(j)));
end
subplot(223);

imagesc(ActValcrop2Qsorted);
xlabel('Time (min)', 'FontSize', fontsize1);
ax = gca;
ax.XTick = 0:400:numtimepoints
ax.XTickLabel = [0:400:numtimepoints] * time_between_frames/60;
ylabel('Worm Number', 'FontSize', fontsize1);
ax = gca;
ax.YTick = 1:numworms2;
%ax.YTickLabel = ActValcrop2totalQsorted* time_between_frames/60;
ax.YTickLabel = group2(index2);  % Comment this out for sequential labeling

set(gca,'FontSize',fontsize2);

fractionwithlessQ1 = [];
fractionwithlessQ2 = [];
maxtime = max([ActValcrop1totalQsorted; ActValcrop2totalQsorted]);
for j=1:maxtime
    fractionwithlessQ1(j) = sum(ActValcrop1totalQsorted < j)/length(group1);
    fractionwithlessQ2(j) = sum(ActValcrop2totalQsorted < j)/length(group2);
end
title('Group 2 (Quiescent = dark)');
% comment next three lines out to remove quiescence labels on right side of heat map
for j=1:numworms2
    text(numtimepoints*1.01,j, num2str(ActValcrop2totalQsorted(j)* time_between_frames/60));
end
%figure(4);clf;
lw=4;
col=[0.2 0.9 0.2];

subplot(222);
plot(fractionwithlessQ1, '-k', 'LineWidth',lw);
xlabel('Time')
ax = gca;
ax.XTick = 0:10:maxtime;
ax.XTickLabel = [0:10:maxtime];
ylabel('CDF ', 'FontSize', fontsize1);
hold on;
plot(fractionwithlessQ2, '-', 'Color',col,'LineWidth',lw);
xlabel('Time (min)', 'FontSize', fontsize1);
ylabel('CDF ', 'FontSize', fontsize1);
ax = gca;
ax.YTick = 0:.2:1;
%title(PathName, 'Interpreter', 'None');
[h,pp] = kstest2(ActValcrop1totalQsorted, ActValcrop2totalQsorted)
title(['Total quiescence, KS Test: p = ' num2str(pp) ' (black=group1)']);
%title('Total quiescence');
set(gca,'FontSize',fontsize2);

%

longestquiescence1 = [];
for j=1:size(ActValcrop1,1)
    d = (ActValcrop1(j,:) >= quiescence_threshold);
    w = [ 1 d 1 ]; % auxiliary vector
    runs_zeros = flip(sort(find(diff(w)==1)-find(diff(w)==-1))); % lengths of runs of 0's
    if length(runs_zeros)>0
        longestquiescence1(j) = runs_zeros(1)* time_between_frames;
    end
end

longestquiescence2 = [];

for j=1:size(ActValcrop2,1)
    d = (ActValcrop2(j,:) >= quiescence_threshold);
    w = [ 1 d 1 ]; % auxiliary vector
    runs_zeros = flip(sort(find(diff(w)==1)-find(diff(w)==-1))); % lengths of runs of 0's
    if length(runs_zeros)>0
        longestquiescence2(j) = runs_zeros(1)* time_between_frames;
    end
end


fractionwithlongestquiescence1shorterthan = [];
fractionwithlongestquiescence2shorterthan = [];
maxtime2 = max([longestquiescence1 longestquiescence2]);
for j=1:maxtime2
    fractionwithlongestquiescence1shorterthan(j) = sum(longestquiescence1 <= j)/length(group1);
    fractionwithlongestquiescence2shorterthan(j) = sum(longestquiescence2 <= j)/length(group2);
end
skip=1;
maxt=skip:skip:maxtime2;
%figure(5);clf;

subplot(224);
plot(maxt,fractionwithlongestquiescence1shorterthan(skip:skip:end), '-k', 'LineWidth',lw);
xlabel('Time (s)', 'FontSize', fontsize1);
ax = gca;
ax.XTick = 0:round(maxtime2/10):maxtime2;
ax.XTickLabel = [0:round(maxtime2/10):maxtime2] * time_between_frames;
ylabel('CDF ', 'FontSize', fontsize1);
hold on;
plot(maxt,fractionwithlongestquiescence2shorterthan(skip:skip:end), '-', 'Color',col, 'LineWidth',lw);
xlabel('Time (s)', 'FontSize', fontsize1);
ylabel('CDF ', 'FontSize', fontsize1);
ax = gca;
ax.YTick = 0:.2:1;
set(gca,'FontSize',fontsize2);
set(gca,'FontSize',fontsize2);
%title(PathName, 'Interpreter', 'None');
[h,p] = kstest2(longestquiescence1, longestquiescence2)
title(['Longest quiescence, KS Test: p = ' num2str(p) ' (black=group1)']);

subplot(221);
text(-2,-3, PathName, 'FontSize', 15);

% Save data in excel format
xlswrite('Group_1',[ActValcrop1totalQsorted longestquiescence1'])
xlswrite('Group_2',[ActValcrop2totalQsorted longestquiescence2'])

%%  CFY ADDED JAN 2019
% 1.  identify transitions from movement to quiescence (a bout of >30 sec). 
% 2. Then assess the length (and speed ?) of the movement bout before the quiescence bout and 
% 3. ask if the duration of the quiescence bout correlated with the duration times the speed of the prior movement bout. 
% 
% Ask next if the correlation is better after training. 

ActValcrop1Qi = [];
durationsQall = [];
durationsNQall = [];

for j=1:numworms1
% for j=1:num
    ActValcrop1Qi = ActValcrop1Q(j,:);
    durations1 = diff(find([1,diff(ActValcrop1Qi),1]));
    if ActValcrop1Qi(1) == 0  % initially quiescent
        durations1 = durations1(2:end);  % crop off initial quiescent period
    end
    if ActValcrop1Qi(end) == 1  % nonquiescent at the end
        durations1 = durations1(1:end-1);  % crop off final nonquiescent period
    end
    durations1Q = durations1(1:2:end);
    durations1NQ = durations1(2:2:end);
    durations1Qall = [durations1Qall durations1Q];
    durations1NQall = [durations1NQall durations1NQ];
end


for j=1:numworms2
    ActValcrop2Qi = ActValcrop2Q(j,:);
    durations2 = diff(find([1,diff(ActValcrop2Qi),1]));
    if ActValcrop2Qi(1) == 0  % initially quiescent
        durations2 = durations2(2:end);  % crop off initial quiescent period
    end
    if ActValcrop2Qi(end) == 1  % nonquiescent at the end
        durations2 = durations2(1:end-1);  % crop off final nonquiescent period
    end
    durations2Q = durations2(1:2:end);
    durations2NQ = durations2(2:2:end);
    durations2Qall = [durations2Qall durations2Q];
    durations2NQall = [durations2NQall durations2NQ];
end


figure(5);
subplot(211);
plot(durations1NQall, durations1Qall, 'o');
xlabel('Active bout length (frames)', 'FontSize', fontsize2);
ylabel('Quiescent bout length (frames)', 'FontSize', fontsize2);
title({[pname(1:end) srcfilename] ; ['Group 1, All bouts included, Correlation R=' num2str(R1(1,2)) ]}, ...
    'Interpreter', 'None');

subplot(212);
plot(durations2NQall, durations2Qall, 'o');
xlabel('Active bout length (frames)', 'FontSize', fontsize2);
ylabel('Quiescent bout length (frames)', 'FontSize', fontsize2);
title(['Group 2, minimum Quiescence=' num2str(durationthreshold) ', Correlation R=' num2str(R2(1,2))]);


durationthreshold = 10;
durations1NQallthresholded = durations1NQall(durations1Qall>durationthreshold);
durations1Qallthresholded = durations1Qall(durations1Qall>durationthreshold);
[R1,P1] = corrcoef(durations1NQallthresholded, durations1Qallthresholded)

durations2NQallthresholded = durations2NQall(durations2Qall>durationthreshold);
durations2Qallthresholded = durations2Qall(durations2Qall>durationthreshold);
[R2,P2] = corrcoef(durations2NQallthresholded, durations2Qallthresholded)

figure(6);
subplot(211);
plot(durations1NQallthresholded, durations1Qallthresholded, 'o');
xlabel('Active bout length (frames)', 'FontSize', fontsize2);
ylabel('Quiescent bout length (frames)', 'FontSize', fontsize2);
title({[pname(1:end) srcfilename] ; ['Group 1, minimum Quiescence=' num2str(durationthreshold) ', Correlation R1=' num2str(R1(1,2))]}, ...
    'Interpreter', 'None');

subplot(212);
plot(durations2NQallthresholded, durations2Qallthresholded, 'o');
xlabel('Active bout length (frames)', 'FontSize', fontsize2);
ylabel('Quiescent bout length (frames)', 'FontSize', fontsize2);
title(['Group 2, minimum Quiescence=' num2str(durationthreshold) ', Correlation R=' num2str(R2(1,2))]);






