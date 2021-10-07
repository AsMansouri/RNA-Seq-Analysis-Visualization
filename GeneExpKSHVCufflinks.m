clc;
clear;

Visibility = 'off';
warning('off');

%% ----Type
% 1- hg19    2- KSHV1 = NC_009333.1(Human)     3- KSHV2 = GQ994935.1(rKSHV)
Type = 1;

IDType = 1; % 1-ENSG  2-ShortID

switch Type
    case 1
        TypeName = 'hg19';
        ResultPath = 'HCC Results/';
        newResults = 'Results/hg19';
        CuffDiffPath = 'CuffDiff/hg19/';
        CufflinkPath = 'Cufflinks_output/hg19/';
    case 2
        TypeName = 'KSHV1';
        ResultPath = 'HCC Results/';
        newResults = 'Results/KSHV1';
        CuffDiffPath = 'CuffDiff/KSHV1/';
        CufflinkPath = 'Cufflinks_output/KSHV1/';
    case 3
        TypeName = 'KSHV2';
        ResultPath = 'HCC Results/';
        newResults = 'Results/KSHV2';
        CuffDiffPath = 'CuffDiff/KSHV2/';
        CufflinkPath = 'Cufflinks_output/KSHV2/';
    otherwise
            disp("-----ERROR----------");
            disp("-----Choose a correct type---------");
end

%% ----GTF annotations 

% Annotobj = GTFAnnotation('D:\Kaposi Sacroma\Data\hg19\GTFs\gencode.v17.annotation.gtf');


%% ----Reading Samples
disp('--------------Reading Samples---------------------');

if exist([TypeName '-Samples-Cufflinks.mat'],'file')
    disp('------------Loading Samples-Cufflinks.mat--------------------');
    load([TypeName '-Samples-Cufflinks.mat'])
else
%     GenesExp = tdfread([ResultPath CuffDiffPath 'genes.fpkm_tracking'],'\t');
    
    Sample(1).Name = 'NC';
%     Sample(2).Name = 'NEG';
    Sample(2).Name = 'GFP';
    Sample(3).Name = 'RGpos';
    for i = 1:length(Sample)
%         if strcmp(Sample(i).Name,'NC') && (Type ~= 1)
%             continue;
%         end
        tmp1 = tdfread([ResultPath Sample(i).Name '/' CufflinkPath 'genes.fpkm_tracking'],'\t');
        Sample(i).GenesIDs = cellstr(tmp1.gene_id);
        %--Order all the Genes and corresponding FPKMs with the First Sample order
        if isempty(Sample(1).GenesIDs)
            [~,~,idxsIntoB] = intersect(Sample(2).GenesIDs,Sample(i).GenesIDs,'stable'); 
        else
            [~,~,idxsIntoB] = intersect(Sample(1).GenesIDs,Sample(i).GenesIDs,'stable');
        end   
        Sample(i).GenesShortIDs = cellstr(tmp1.gene_short_name(idxsIntoB,:));
        Sample(i).GenesFPKMs = tmp1.FPKM(idxsIntoB);
        Sample(i).GenesIDs = cellstr(tmp1.gene_id(idxsIntoB,:));
        Sample(i).GenesENSG = extractBefore(Sample(i).GenesIDs,'.');
    end
    
    disp('------------Saving Samples-Cufflinks.mat--------------------');
    save([TypeName '-Samples-Cufflinks.mat'],'Sample')
end

%ImportantVars = {'Visibility','Sample','ImportantVars'};
%clearvars('-except',ImportantVars{:})
disp('--------------------------------------------------');

if ~exist(newResults,'dir')
    mkdir(newResults);
end
newResults = [newResults '/'];


disp('----------Saving each Samples FPKMs Values------------');
% for i = 1:length(Sample)
%     %-----FPKMs Values
%     figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%     bar(Sample(i).GenesFPKMs);
%     title(Sample(i).Name);
%     if Type ~= 1
%         set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
%     end
%     xlabel('Genes')
%     ylabel('FPKM')
%     saveas(gcf,[newResults Sample(i).Name '-FPKMs.jpg']);
%     saveas(gcf,[newResults Sample(i).Name '-FPKMs.fig']);
%     SavePDF(gcf,[newResults Sample(i).Name '-FPKMs']);
%     close all force;
%     %----log2 of FPKMs Values
%     figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%     bar(log2(Sample(i).GenesFPKMs+1));
%     title(Sample(i).Name);
%     if Type ~= 1
%         set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
%     end
%     xlabel('Genes')
%     ylabel('Log2 FPKM')
%     saveas(gcf,[newResults Sample(i).Name '-Log2 FPKMs.jpg']);
%     saveas(gcf,[newResults Sample(i).Name '-Log2 FPKMs.fig']);
%     SavePDF(gcf,[newResults Sample(i).Name '-Log2 FPKMs']);
%     close all force;
% end

%-----FPKMs Values of all samples in one figure
% figure('Visible',Visibility,'Position', get(0, 'Screensize'))
% for i = 1:length(Sample)
%     if i == 1
%         title('Gene Coverage');
%     end
%     subplot(length(Sample),1,i)
%     bar(Sample(i).GenesFPKMs);
% %     xlabel('Genes')
%     ylabel(Sample(i).Name)
% end
% if Type ~= 1
%     set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
% end
% saveas(gcf,[newResults 'All -FPKMs.jpg']);
% saveas(gcf,[newResults 'All -FPKMs.fig']);
% SavePDF(gcf,[newResults 'All -FPKMs']);
% close all force;

% figure('Visible',Visibility,'Position', get(0, 'Screensize'))
% for i = 1:length(Sample)
%     if i == 1
%         title('Gene Coverage (log2)');
%     end
%     subplot(length(Sample),1,i)
%     bar(log2(Sample(i).GenesFPKMs+1));
% %     xlabel('Genes')
%     ylabel(Sample(i).Name)
% end
% if Type ~= 1
%     set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
% end
% saveas(gcf,[newResults 'All -Log2 FPKMs.jpg']);
% saveas(gcf,[newResults 'All -Log2 FPKMs.fig']);
% SavePDF(gcf,[newResults 'All -Log2 FPKMs']);
% close all force;

disp('--------------------------------------------------');

%% ----- Filtering UnExpressed genes at all groups

%remove noises
FPKM_Thresh = 1;

CommGeneFPKMsMat = [];
CommGeneName = {};

for i = 1:length(Sample)
    CommGeneFPKMsMat =  [CommGeneFPKMsMat; Sample(i).GenesFPKMs'];
    CommGeneName(end+1) = {Sample(i).Name};
end

CommGeneFPKMsMat(CommGeneFPKMsMat < FPKM_Thresh) = 0;

tmp = find(sum(CommGeneFPKMsMat)==0);

for i = 1:length(Sample)
    Sample(i).GenesIDs(tmp) = [];
    Sample(i).GenesENSG(tmp) = [];
    Sample(i).GenesFPKMs(tmp) = [];
    Sample(i).GenesShortIDs(tmp) = [];
end
CommGeneFPKMsMat(:,tmp) = [];

%% ----- Bars Plots All Smples
% figure('Visible',Visibility,'Position', get(0, 'Screensize'))
% bar(CommGeneFPKMsMat')
% title('All Strains');
% legend(CommGeneName)
% if Type ~= 1
%     set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
% end
% xlabel('Genes')
% ylabel('FPKMs')
% saveas(gcf,[newResults 'All.jpg']);
% saveas(gcf,[newResults 'All.fig']);
% SavePDF(gcf,[newResults 'All']);
% close all force;
%  
% figure('Visible',Visibility,'Position', get(0, 'Screensize'))
% bar(log2(CommGeneFPKMsMat'))
% title('All Strains');
% legend(CommGeneName)
% if Type ~= 1
%     set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
% end
% xlabel('Genes')
% ylabel('log2 FPKMs')
% saveas(gcf,[newResults 'All - log2.jpg']);
% saveas(gcf,[newResults 'All - log2.fig']);
% SavePDF(gcf,[newResults 'All - log2']);
% close all force;


%% ---- Pairwise regulation for KSHV & HEATMAP

% if Type ~= 1
%     Folds = {};
%     for i = 2:length(Sample)
%         for j = i+1:length(Sample)
%             noise = 0.0001;
%             tmp1FPKMs = Sample(i).GenesFPKMs;
%             tmp2FPKMs = Sample(j).GenesFPKMs;
%             Folds(end+1).Name = [Sample(i).Name ' / ' Sample(j).Name];
%             Folds(end).Values = (tmp1FPKMs+noise) ./ (tmp2FPKMs+noise);
%             ind = find(Folds(end).Values > 2 | Folds(end).Values < 0.5);
%             figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%             bar(log2(Folds(end).Values(ind)))
%             title({Folds(end).Name ;'Folds change'});
%             if Type ~= 1
%                 set(gca,'XTick',[1:length(Sample(i).GenesShortIDs(ind))],'XTickLabel',Sample(i).GenesShortIDs(ind),'XTickLabelRotation',90);
%             end
%             xlabel('Genes')
%             ylabel('log2 (Fold Change)')
%             saveas(gcf,[newResults Sample(i).Name ' - ' Sample(j).Name ' Folds change - log2.jpg']);
%             saveas(gcf,[newResults Sample(i).Name ' - ' Sample(j).Name ' Folds change - log2.fig']);
%             SavePDF(gcf,[newResults Sample(i).Name ' - ' Sample(j).Name ' Folds change - log2']);
%             
%             close all force;
% 
%             tmp1FPKMs = Sample(j).GenesFPKMs;
%             tmp2FPKMs = Sample(i).GenesFPKMs;
%             Folds(end+1).Name = [Sample(j).Name ' / ' Sample(i).Name];
%             Folds(end).Values = (tmp1FPKMs+noise) ./ (tmp2FPKMs+noise);
%             ind = find(Folds(end).Values > 2 | Folds(end).Values < 0.5);
%             figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%             bar(log2(Folds(end).Values(ind)))
%             title({Folds(end).Name ;'Folds change'});
%             if Type ~= 1
%                 set(gca,'XTick',[1:length(Sample(i).GenesShortIDs(ind))],'XTickLabel',Sample(i).GenesShortIDs(ind),'XTickLabelRotation',90);
%             end
%             xlabel('Genes')
%             ylabel('log2 (Fold Change)')
%             saveas(gcf,[newResults Sample(j).Name ' - ' Sample(i).Name ' Folds change - log2.jpg']);
%             saveas(gcf,[newResults Sample(j).Name ' - ' Sample(i).Name ' Folds change - log2.fig']);
%             SavePDF(gcf,[newResults Sample(j).Name ' - ' Sample(i).Name ' Folds change - log2']);
%             close all force;
%         end
%     end
%     
%     c=zeros(length(Sample),3);
%     c([1], :) = ones(1,3).*[0,0,1]; % NC
%     c([2], :) = ones(1,3).*[0,0,0]; % NEG
%     c([3], :) = ones(1,3).*[0,1,0]; % GFP
%     c([4], :) = ones(1,3).*[1,0,0]; % RGpos
% 
%     for i = 2:length(Sample)
%         tmpFolds = [];
%         Samps = [2:4];
%         Samps(Samps == i) = [];
%         tmpNames = {Sample(i).Name};
%         for j = 1:length(Samps)
%             noise = 0.0001;
%             tmp1FPKMs = Sample(Samps(j)).GenesFPKMs;
%             tmp2FPKMs = Sample(i).GenesFPKMs;
%             tmpNames(end+1) = {Sample(Samps(j)).Name};
%             tmpFolds(:,j) = (tmp1FPKMs+noise) ./ (tmp2FPKMs+noise); 
%         end
%         figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%         hb = bar(log2(tmpFolds));
%         for c1 = 1:length(Samps)
%             set(hb(c1), 'FaceColor',c(Samps(c1),:))
%         end
%         legend(tmpNames(2:end))
%         title({['( ' tmpNames{2} ' & ' tmpNames{3} ') / ' tmpNames{1}] ;'Folds change'});
%         if Type ~= 1
%             set(gca,'XTick',[1:length(Sample(i).GenesShortIDs)],'XTickLabel',Sample(i).GenesShortIDs,'XTickLabelRotation',90);
%         end
%         xlabel('Genes')
%         ylabel('log2 (Fold Change)')
%         saveas(gcf,[newResults '( ' tmpNames{2} ' & ' tmpNames{3} ') - ' tmpNames{1} ' Folds change - log2.jpg']);
%         saveas(gcf,[newResults '( ' tmpNames{2} ' & ' tmpNames{3} ') - ' tmpNames{1} ' Folds change - log2.fig']);
%         SavePDF(gcf,[newResults '( ' tmpNames{2} ' & ' tmpNames{3} ') - ' tmpNames{1} ' Folds change - log2']);
%         close all force;
%     end
%     
%     figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%     h = heatmap(Sample(i).GenesShortIDs,{'NEG','GFP','RGpos'},CommGeneFPKMsMat(2:end,:));
%     h.ColorScaling = 'scaledcolumns';
%     h.Title = 'Gene expressions';
%     h.YLabel = 'Samples';
%     h.XLabel = 'Genes';
%     h.Colormap = parula;
%     saveas(gcf,[newResults 'Heatmap all fpkms.jpg']);
%     saveas(gcf,[newResults 'Heatmap all fpkms.fig']);
%     SavePDF(gcf,[newResults 'Heatmap all fpkms']);
%     close all force;
%     
%     figure('Visible',Visibility,'Position', get(0, 'Screensize'))
%     h = heatmap(Sample(i).GenesShortIDs,{'NEG','GFP','RGpos'},log2(CommGeneFPKMsMat(2:end,:)+1));
%     h.Title = 'Log2 of Gene expressions';
%     h.YLabel = 'Samples';
%     h.XLabel = 'Genes';
%     h.Colormap = parula;
%     saveas(gcf,[newResults 'Heatmap all fpkms- log2.jpg']);
%     saveas(gcf,[newResults 'Heatmap all fpkms- log2.fig']);
%     SavePDF(gcf,[newResults 'Heatmap all fpkms- log2']);
%     close all force;
%     
% end

% ---- Genes Groups

FPKM_Thresh = 1;
NoiseRatio = 0.8; % ratio of accepted replicates

if Type == 1
    if IDType == 1
        [Genes] = GroupingGenes (Sample(1).GenesENSG, CommGeneFPKMsMat, FPKM_Thresh, NoiseRatio);
    elseif IDType == 2
        [Genes] = GroupingGenes (Sample(1).GenesShortIDs, CommGeneFPKMsMat, FPKM_Thresh, NoiseRatio);
    end
else
    [Genes] = GroupingGenes (Sample(1).GenesShortIDs, CommGeneFPKMsMat, FPKM_Thresh, NoiseRatio);
end



%% ----- BoxPlot
BoxPlots (CommGeneFPKMsMat, CommGeneName, newResults, Visibility);


%% ----- SVD

% SVDs(CommGeneFPKMsMat, CommGeneName, '2D', newResults, Visibility)
% SVDs(CommGeneFPKMsMat, CommGeneName, '3D', newResults, Visibility)


%% ----- Cluster gram
try 
    ClusterGrams (CommGeneFPKMsMat, CommGeneName, Sample(1).GenesShortIDs, Type, newResults, ' ');

    if Type ~= 1
        ClusterGrams (CommGeneFPKMsMat(2:end,:), CommGeneName(2:end), Sample(1).GenesShortIDs, Type, newResults, ' Infected Samples');
    end
catch
    disp('No Clustergram');
end
%% ---- venn diagram

% if Type == 1
%     figure, axis equal, axis off
%     A = [length(Genes.RGpos) length(Genes.GFP) length(Genes.NEG)];
%     I = [length(Genes.GFP_RGpos) length(Genes.NEG_RGpos) length(Genes.NEG_GFP) ...
%          length(Genes.RGpos_GFP_NEG)];
%     [~,S] = venn(A,I,'FaceColor',{'r','g','k'},'EdgeColor','black','ErrMinMode', 'ChowRodgers');
%     for i = 1:length(S.ZoneCentroid)
%         text(S.ZoneCentroid(i,1),S.ZoneCentroid(i,2),num2str(S.ZonePop(i)));
%     end
%     title({'Venn Diagram';'RGpos & GFP & NEG'})
%     saveas(gcf,[newResults 'Venn - RGpos & GFP & NEG.jpg']);
%     saveas(gcf,[newResults 'Venn - RGpos & GFP & NEG.jpg']);
%     close all force;
% else
%     disp('-----  No Venn Diagram -----');
% end



%% ---Save genes name for each group in an EXCEL file

GroupsGenesEXCEL (Genes, 'Genes', 'Genes', newResults, Visibility);

AllGenesExpressionEXCEL (Sample, CommGeneFPKMsMat, 'All Genes Expression', newResults);

%% EachGroup Expressed Genes
SamplesIndx = [1 ; 2 ; 3];

%--Just NC
tmpL = 'Just NC';
if Type == 1
    if IDType == 1
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesENSG,Genes.Just_NC,'stable');
    elseif IDType == 2
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_NC,'stable');
    end
else
    [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_NC,'stable');
end
tmpGeneFPKMs = CommGeneFPKMsMat(SamplesIndx(1),RefindxGenes);

GroupSpecExpGenes(tmpGeneFPKMs, GeneENSG, tmpL, Type, newResults, Visibility);



%--Just GFP
tmpL = 'Just GFP';
if Type == 1
    if IDType == 1
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesENSG,Genes.Just_GFP,'stable');
    elseif IDType == 2
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_GFP,'stable');
    end
else
    [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_GFP,'stable');
end
tmpGeneFPKMs = CommGeneFPKMsMat(SamplesIndx(2),RefindxGenes);

GroupSpecExpGenes(tmpGeneFPKMs, GeneENSG, tmpL, Type, newResults, Visibility);


%--Just RGpos
tmpL = 'Just RGpos';
if Type == 1
    if IDType == 1
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesENSG,Genes.Just_RGpos,'stable');
    elseif IDType == 2
        [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_RGpos,'stable');
    end
else
    [GeneENSG,RefindxGenes,~] = intersect(Sample(1).GenesShortIDs,Genes.Just_RGpos,'stable');
end
tmpGeneFPKMs = CommGeneFPKMsMat(SamplesIndx(3),RefindxGenes);

GroupSpecExpGenes(tmpGeneFPKMs, GeneENSG, tmpL, Type, newResults, Visibility);


% TopNexpressed(Sample,[2,3], 50, 'TopExpressed.xlsx', 'Top ', newResults, Visibility)

%% -------CuffDiff Analysis
disp('------------Loading CuffDiff "gene_exp.diff"--------------------');

if exist([TypeName ' - CuffDiffTable.mat'],'file')
    disp('------------Loading CuffDiffTable.mat--------------------');
    load([TypeName ' - CuffDiffTable.mat'])
    disp('--------------------------------------------------');
else
    disp('------------Saving CuffDiffTable.mat--------------------');
    GeneExpCuffDiff = tdfread([ResultPath CuffDiffPath 'gene_exp.diff'],'\t');

    tmp=strtrim(cellstr(GeneExpCuffDiff.gene_id));

    index = find(ismember(tmp, Sample(1).GenesIDs)); % Find expressed genes 


    T1 = table(strtrim(GeneExpCuffDiff.gene_id(index,:)),strtrim(GeneExpCuffDiff.gene(index,:)),strtrim(GeneExpCuffDiff.locus(index,:)),...
           strtrim(GeneExpCuffDiff.sample_1(index,:)),strtrim(GeneExpCuffDiff.sample_2(index,:)),GeneExpCuffDiff.value_1(index,:),...
           GeneExpCuffDiff.value_2(index,:),GeneExpCuffDiff.log20x28fold_change0x29(index,:),GeneExpCuffDiff.p_value(index,:),...
           GeneExpCuffDiff.q_value(index,:),GeneExpCuffDiff.significant(index,:),...
        'VariableNames',{'ENSAMBLE','Gene','Locus','Sample_1','Sample_2','Value_1','Value_2','Fold','p_value','q_value','Significant'});

    save([TypeName ' - CuffDiffTable.mat'],'T1')
    disp('--------------------------------------------------');
end

%--Check the p-value q-value
tmp = table2array(T1(:,11));
T1(tmp(:,1)=='n' ,:) = [];
tmp = table2array(T1(:,4));

if Type==1
    T1(tmp(:,1)~='N',:) = [];
%     T1(tmp(:,2)~='C',:) = [];
end

%-Fold Table is ready. NC vs others

%--Check the fold change > 2
tmp = table2array(T1(:,8));
T1(tmp >= 1 & tmp <= -1 ,:) = [];



CLabels = cellstr(strtrim(unique(table2cell(T1(:,4)))));
SLabels = cellstr(strtrim(unique(table2cell(T1(:,5)))));
% [SLabels(1), SLabels(2)] = deal(SLabels(2), SLabels(1)); % Swap date1 and date2
SLabels = sort(SLabels);

Labels = [CLabels;SLabels];    %% => CHECK their order:  NC,GFP,RGpos 



if IDType == 1
    tmpUniqgenes = extractBefore(cellstr(strtrim(unique(table2cell(T1(:,1))))),'.');
    tmpgenes = extractBefore(cellstr(strtrim(table2cell(T1(:,1)))),'.');
else
    tmpUniqgenes = cellstr(strtrim(unique(table2cell(T1(:,2)))));
    tmpgenes = cellstr(strtrim(table2cell(T1(:,2))));
end
    
    


FoldsValuesMat = zeros(length(tmpUniqgenes),length(CLabels)*length(SLabels));
FPKMsValuesMat = zeros(length(tmpUniqgenes),length(CLabels)+length(SLabels));
FDRValuesMat = zeros(length(tmpUniqgenes),length(CLabels)*length(SLabels));
FoldsValuesMat = ones(length(tmpUniqgenes),length(CLabels)*length(SLabels))*(-100);
FPKMsValuesMat = ones(length(tmpUniqgenes),length(CLabels)+length(SLabels))*(-1);
FDRValuesMat = ones(length(tmpUniqgenes),length(CLabels)*length(SLabels));
% 
S1tmp = cellstr(table2cell(T1(:,4)));
S2tmp = cellstr(table2cell(T1(:,5)));
S1FPKMstmp = table2array(T1(:,6));
S2FPKMstmp = table2array(T1(:,7));
S12Foldstmp = table2array(T1(:,8));
FDRstmp = table2array(T1(:,10));

% Table_Under_Exp_Combined = table();

for i = 1:length(tmpUniqgenes)
    index = find(ismember(tmpgenes, tmpUniqgenes(i)));
    for ind = 1:length(index)
        indS1 = find(ismember(CLabels, strtrim(S1tmp(index(ind)))));
        indS2 = find(ismember(SLabels, strtrim(S2tmp(index(ind)))));
        FPKMsValuesMat(i,indS1)=S1FPKMstmp(index(ind));
        FPKMsValuesMat(i,indS2+length(CLabels))=S2FPKMstmp(index(ind));
        FoldsValuesMat(i,((indS2-1)*length(CLabels)+indS1))=S12Foldstmp(index(ind));
        FDRValuesMat(i,((indS2-1)*length(CLabels)+indS1))=FDRstmp(index(ind));
    end
end

if IDType == 2
    [~,~,RefB] = intersect(tmpUniqgenes,Sample(1).GenesShortIDs,'stable');
    tmpUniqgenesShortID = Sample(1).GenesShortIDs(RefB);
else
    [~,~,RefB] = intersect(tmpUniqgenes,extractBefore(Sample(1).GenesIDs,'.'),'stable');
    tmpUniqgenesShortID = Sample(1).GenesShortIDs(RefB);
end

T = table(tmpUniqgenes,tmpUniqgenesShortID,...
    FPKMsValuesMat(:,1),FPKMsValuesMat(:,2),FPKMsValuesMat(:,3),...
    FoldsValuesMat(:,1),FoldsValuesMat(:,2),...
    FDRValuesMat(:,1),FDRValuesMat(:,2));

xlswrite([newResults 'SignifRegulatedTable.xlsx'],table2cell(T));

tmpFPKMsMat=sortrows(FPKMsValuesMat,[1:size(FPKMsValuesMat,2)]);

    
for s = 1:length(SLabels)
    tmpindex = find(ismember(strtrim(S2tmp),SLabels(s)));
    RegulatedSamples(s).Name = SLabels(s);
    RegulatedSamples(s).GeneIDs = tmpgenes(tmpindex);
    RegulatedSamples(s).FPKMsS1 = S1FPKMstmp(tmpindex);
    RegulatedSamples(s).FPKMsS2 = S2FPKMstmp(tmpindex);
    RegulatedSamples(s).Folds = S12Foldstmp(tmpindex);
    RegulatedSamples(s).FPKMSFoldsMat = [S1FPKMstmp(tmpindex) S2FPKMstmp(tmpindex) S12Foldstmp(tmpindex)];
end

%------All type of regulation:
RegulatedGenes.GFP = RegulatedSamples(1).GeneIDs;
RegulatedGenes.RGpos = RegulatedSamples(2).GeneIDs;

RegulatedGenes.GFP_RGpos = RegulatedGenes.GFP(ismember(RegulatedGenes.GFP,RegulatedGenes.RGpos));

RegulatedGenes.Just_GFP = RegulatedGenes.GFP(~ismember(RegulatedGenes.GFP,RegulatedGenes.RGpos));
RegulatedGenes.Just_RGpos = RegulatedGenes.RGpos(~ismember(RegulatedGenes.RGpos,RegulatedGenes.GFP));

[RegulatedFoldsValuesMat] = RegulatedFolds (RegulatedGenes.GFP_RGpos, RegulatedSamples, [1,2]);
[RegulatedGenes] = RegulationBars(RegulatedGenes, RegulatedFoldsValuesMat, RegulatedGenes.GFP_RGpos, ['g','r'], 'GFP_RGpos', {'GFP','RGpos'}, 'GFP & RGpos Regulations- ', newResults, Visibility);

[RegulatedFoldsValuesMat] = RegulatedFolds (RegulatedGenes.Just_GFP, RegulatedSamples, [1]);
[RegulatedGenes] = RegulationBars(RegulatedGenes, RegulatedFoldsValuesMat, RegulatedGenes.Just_GFP, ['g'], 'Just_GFP', {'GFP'}, 'Just GFP Regulations- ', newResults, Visibility);

[RegulatedFoldsValuesMat] = RegulatedFolds (RegulatedGenes.Just_RGpos, RegulatedSamples, [2]);
[RegulatedGenes] = RegulationBars(RegulatedGenes, RegulatedFoldsValuesMat, RegulatedGenes.Just_RGpos, ['r'], 'Just_RGpos', {'RGpos'}, 'Just RGpos Regulations- ', newResults, Visibility);

%% ---Save Regulated genes name for each group in an EXCEL file

GroupsGenesEXCEL (RegulatedGenes, 'RegulatedGenes0', 'RegulatedGenes0', newResults, Visibility);


% Genes = RegulatedGenes;
% name = 'RegulatedGenes';
% Titles = 'RegulatedGenes';
% Path = newResults;

%%

%% C) GSE84237 - Sychev,Laguoff
OthersSigSamples={};

[~,~,raw1] = xlsread('GSEs\GSE84237_TIME_KSHVvsTIME_Mock-Trimmed.xlsx');

OthersSigSamples(1).Name = 'GSE84237';
OthersSigSamples(end).GeneIDs = raw1(2:end,1);
OthersSigSamples(end).GeneShortIDs = raw1(2:end,2);
OthersSigSamples(end).Folds = [raw1{2:end,3}];
OthersSigSamples(end).FDR = [raw1{2:end,6}];
OthersSigSamples(end).signif = [raw1{2:end,7}];

tmpindx = OthersSigSamples(1).FDR < 0.01;
OthersSigSamples(2).Name = 'GSE84237-Signif-FDR';
OthersSigSamples(end).GeneIDs = OthersSigSamples(1).GeneIDs(tmpindx);
OthersSigSamples(end).GeneShortIDs = OthersSigSamples(1).GeneShortIDs(tmpindx);
OthersSigSamples(end).Folds = OthersSigSamples(1).Folds(tmpindx);
OthersSigSamples(end).FDR = OthersSigSamples(1).FDR(tmpindx);
OthersSigSamples(end).signif = OthersSigSamples(1).signif(tmpindx);


tmpindx = OthersSigSamples(1).signif > 0 ;
OthersSigSamples(3).Name = 'GSE84237-Signif-sig';
OthersSigSamples(end).GeneIDs = OthersSigSamples(1).GeneIDs(tmpindx);
OthersSigSamples(end).GeneShortIDs = OthersSigSamples(1).GeneShortIDs(tmpindx);
OthersSigSamples(end).Folds = OthersSigSamples(1).Folds(tmpindx);
OthersSigSamples(end).FDR = OthersSigSamples(1).FDR(tmpindx);
OthersSigSamples(end).signif = OthersSigSamples(1).signif(tmpindx);


tmpindx = OthersSigSamples(3).Folds > 0 ;
OthersSigSamples(4).Name = 'GSE84237-Signif-sig-Up';
OthersSigSamples(end).GeneIDs = OthersSigSamples(3).GeneIDs(tmpindx);
OthersSigSamples(end).GeneShortIDs = OthersSigSamples(3).GeneShortIDs(tmpindx);
OthersSigSamples(end).Folds = OthersSigSamples(3).Folds(tmpindx);
OthersSigSamples(end).FDR = OthersSigSamples(3).FDR(tmpindx);
OthersSigSamples(end).signif = OthersSigSamples(3).signif(tmpindx);

tmpindx = OthersSigSamples(3).Folds < 0 ;
OthersSigSamples(5).Name = 'GSE84237-Signif-sig-Down';
OthersSigSamples(end).GeneIDs = OthersSigSamples(3).GeneIDs(tmpindx);
OthersSigSamples(end).GeneShortIDs = OthersSigSamples(3).GeneShortIDs(tmpindx);
OthersSigSamples(end).Folds = OthersSigSamples(3).Folds(tmpindx);
OthersSigSamples(end).FDR = OthersSigSamples(3).FDR(tmpindx);
OthersSigSamples(end).signif = OthersSigSamples(3).signif(tmpindx);


IntersectGFP_RG_TIME = {};
IntersectGFP_RG_TIME(1).Name = 'GFPRG_Up_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.GFP_RGpos_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(2).Name = 'GFPRG_Up_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.GFP_RGpos_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(3).Name = 'GFPRG_Down_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.GFP_RGpos_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(4).Name = 'GFPRG_Down_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.GFP_RGpos_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(5).Name = 'GFP_Down_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.Just_GFP_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];


IntersectGFP_RG_TIME(6).Name = 'GFP_Down_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.Just_GFP_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(7).Name = 'GFP_Up_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.Just_GFP_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(8).Name = 'GFP_Up_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.Just_GFP_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(9).Name = 'RG_Down_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.Just_RGpos_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(10).Name = 'RG_Down_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.Just_RGpos_DownReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(11).Name = 'RG_Up_TIME_Up';
[~,tmpindx,~] = intersect(OthersSigSamples(4).GeneIDs,RegulatedGenes.Just_RGpos_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(4).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(4).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(4).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(4).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(12).Name = 'RG_Up_TIME_Down';
[~,tmpindx,~] = intersect(OthersSigSamples(5).GeneIDs,RegulatedGenes.Just_RGpos_UpReg,'stable');
IntersectGFP_RG_TIME(end).GeneIDs = OthersSigSamples(5).GeneIDs(tmpindx);
IntersectGFP_RG_TIME(end).GeneShortIDs = OthersSigSamples(5).GeneShortIDs(tmpindx);
IntersectGFP_RG_TIME(end).Folds = OthersSigSamples(5).Folds(tmpindx)';
IntersectGFP_RG_TIME(end).FDR = OthersSigSamples(5).FDR(tmpindx);
[~,~,tmpindx] = intersect(IntersectGFP_RG_TIME(end).GeneIDs,tmpUniqgenes,'stable');
IntersectGFP_RG_TIME(end).AllFolds = [FoldsValuesMat(tmpindx,:),IntersectGFP_RG_TIME(end).Folds];

IntersectGFP_RG_TIME(13).Name = 'AllCombined';
tmpgenes = {};
tmpshortID = {};
tmpAllFolds = [];
for i = 1:length(IntersectGFP_RG_TIME)-1
    tmpgenes = [tmpgenes;IntersectGFP_RG_TIME(i).GeneIDs];
    tmpshortID = [tmpshortID;IntersectGFP_RG_TIME(i).GeneShortIDs];
    tmpAllFolds = [tmpAllFolds;IntersectGFP_RG_TIME(i).AllFolds];
end
IntersectGFP_RG_TIME(end).GeneIDs = tmpgenes;
IntersectGFP_RG_TIME(end).GeneShortIDs = tmpshortID;
IntersectGFP_RG_TIME(end).AllFolds = tmpAllFolds;


%%
[ia,ib,ic] = intersect(OthersSigSamples(3).GeneIDs,RegulatedSamples(2).GeneIDs,'stable');
% [a,b,c] = intersect(OthersSigSamples(3).GeneIDs,RegulatedSamples(2).GeneIDs,'stable');
% [a,b,c] = intersect(OthersSigSamples(3).GeneIDs,RegulatedSamples(3).GeneIDs,'stable');

%% ?

% if IDType == 1
%     tmpUniqgenes = extractBefore(cellstr(strtrim(unique(table2cell(T1(:,1))))),'.');
%     tmpgenes = extractBefore(cellstr(strtrim(table2cell(T1(:,1)))),'.');
% else
%     tmpUniqgenes = cellstr(strtrim(unique(table2cell(T1(:,2)))));
%     tmpgenes = cellstr(strtrim(table2cell(T1(:,2))));
% end
%     
% FoldsValuesMatAll = ones(length(tmpUniqgenes),length(CLabels)*length(SLabels))*(-20);
% FPKMsValuesMatAll = zeros(length(tmpUniqgenes),length(CLabels)+length(SLabels));
% 
% S1tmp = cellstr(table2cell(T1(:,4)));
% S2tmp = cellstr(table2cell(T1(:,5)));
% S1FPKMstmp = table2array(T1(:,6));
% S2FPKMstmp = table2array(T1(:,7));
% S12Foldstmp = table2array(T1(:,8));
% 
% for i = 1:length(tmpUniqgenes)
%     index = find(ismember(tmpgenes, tmpUniqgenes(i)));
%     for ind = 1:length(index)
%         indS1 = find(ismember(CLabels, strtrim(S1tmp(index(ind)))));
%         indS2 = find(ismember(SLabels, strtrim(S2tmp(index(ind)))));
%         FPKMsValuesMatAll(i,indS1)=S1FPKMstmp(index(ind));
%         FPKMsValuesMatAll(i,indS2+length(CLabels))=S2FPKMstmp(index(ind));
%         FoldsValuesMatAll(i,((indS2-1)*length(CLabels)+indS1))=S12Foldstmp(index(ind));
%     end
% end

%% show the clustergram of Folds,

%filter inf and -inf
% FoldsValuesMatAll(FoldsValuesMatAll > 20) = 20;
% FoldsValuesMatAll(FoldsValuesMatAll < -20) = -20;
% 
% ClusterGrams (FoldsValuesMatAll', SLabels, Sample(1).GenesShortIDs, Type, newResults, ' Significant Differentially Exp');



%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%Pathways%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
%% ---Save genes name for each group in an EXCEL file

if exist('GeneIDConvertor.mat','file')
    disp(['------------Loading ' 'GeneIDConvertor.mat' '  Pathways info  --------------------']);
    load('GeneIDConvertor.mat')
else
    [~,~,raw1] = xlsread('GenesConvertor.xlsx','NCBI to ENSG');
    [~,~,raw2] = xlsread('GenesConvertor.xlsx','KEGG 2 NCBI');
    [a,b,c] = intersect(raw1(2:end,1),raw2(2:end,1),'stable');
    raw1(b+1,4) = raw2(c+1,2);
    raw1(1,4) = {'KEGG'};
    raw1(1,1) = {'RefSeq_Protein'};
    raw1(1,2) = {'ENSG'};
    GeneIDConverter = raw1(:,[1 2 4 3]);
    tmp = {};
    for i = 2:size(GeneIDConverter,1)
        tmp{i,1} = num2str(GeneIDConverter{i,3});
    end
    GeneIDConverter{2:end,3} = tmp{2:end,1};
    save('GeneIDConvertor.mat','GeneIDConverter');
end

%% ---- Pathways

mkdir('Results\Pathways');

if exist('Results\Pathways\AllPathwaysInfo.mat','file')
    disp(['------------Loading ' 'Results\Pathways\AllPathwaysInfo.mat' '  Pathways info  --------------------']);
    load('Results\Pathways\AllPathwaysInfo.mat')
else
    Org = 'hsa';
 
    Pathwaysstr = strsplit(webread(['http://rest.kegg.jp/list/pathway/' Org]),'\n');
    Pathways = cell(length(Pathwaysstr)-1,1);

    AllMainClasses = {};
    AllClasses = {};
    for p = 1:length(Pathwaysstr)-1
        tmp = strsplit(Pathwaysstr{p},'\t');
        Pathways{p}.Entry = extractAfter(tmp{1},':');
        Pathways{p}.Name = extractBefore(tmp{2},' - ');
        Pathways{p}.Organism = extractAfter(tmp{2},' - ');
        tmp = strsplit(webread(['http://rest.kegg.jp/get/' Pathways{p}.Entry]),'\n');
        tmpIndx = find(contains(tmp,'CLASS')); 
        if ~isempty(tmpIndx)
            Pathways{p}.MainClass = strtrim(extractBefore(strtrim(extractAfter(tmp{tmpIndx},' ')),';'));
            if ~contains(Pathways{p}.MainClass,AllMainClasses)
                AllMainClasses{end+1} = Pathways{p}.MainClass;
            end
            Pathways{p}.Class = strtrim(extractAfter(strtrim(extractAfter(tmp{tmpIndx},' ')),';'));
            if ~contains(Pathways{p}.Class,AllClasses)
                AllClasses{end+1} = Pathways{p}.Class;
            end
        else
            Pathways{p}.MainClass = 'Not Defined';
            Pathways{p}.Class = 'Not Defined';
        end
        tmpIndx = find(contains(tmp,'GENE'));
        if ~isempty(tmpIndx)
            PathGenes = {};
            tmp{tmpIndx}(1:4) = [];
            while tmp{tmpIndx}(1) == ' '
                tmpg = strtrim(tmp{tmpIndx});
                PathGenes(end+1).Entry = extractBefore(strtrim(tmpg),' ');
                PathGenes(end).Name = extractAfter(extractBefore(strtrim(tmpg),';'),' '); 
                PathGenes(end).Describtion = extractAfter(strtrim(tmpg),' ');
                tmpIndx = tmpIndx + 1; 
            end
            Pathways{p}.Genes = PathGenes;
        else
            Pathways{p}.Genes = 'No Genes';
        end

    end

    ConvGenes = strsplit(webread(['http://rest.kegg.jp/conv/' Org '/ncbi-proteinid']),'\n')';
    GeneConvertor = {};
    for p = 1:length(ConvGenes)-1
        tmp = strsplit(ConvGenes{p},'\t');
        GeneConvertor.NCBI{p} = extractAfter(tmp{1},':');
        GeneConvertor.KEGG{p} = extractAfter(tmp{2},':');    
    end

    tmpp = [];
    AllPathwaysInfo = {};
    for p = 1:length(Pathways)
        if ~isstruct(Pathways{p}.Genes)
            tmpp = [tmpp p];
        else
            AllPathwaysInfo{end+1} = Pathways{p};

        end
    end
%     ConvGenes = strsplit(webread(['http://rest.kegg.jp/conv/' Org '/ncbi-proteinid']),'\n')';
%     GeneConvertor1 = {};
%     for p = 1:length(ConvGenes)-1
%         tmp = strsplit(ConvGenes{p},'\t');
%         GeneConvertor1.NCBI{p} = extractAfter(tmp{1},':');
%         GeneConvertor1.KEGG{p} = extractAfter(tmp{2},':');    
%     end
    
    for p = 1:length(AllPathwaysInfo)
        [a,~,c] = intersect({AllPathwaysInfo{p}.Genes.Entry},GeneIDConverter(:,3),'stable');
        AllPathwaysInfo{p}.GenesShortID = {AllPathwaysInfo{p}.Genes.Name};
        AllPathwaysInfo{p}.GenesENSG = GeneIDConverter(c,2);
        AllPathwaysInfo{p}.GenesRefSeq = GeneIDConverter(c,1);
        AllPathwaysInfo{p}.GenesName = GeneIDConverter(c,4);
    end    
    
    for p = 1:length(AllPathwaysInfo)
        [a,b,c] = intersect(Sample(1).GenesENSG,AllPathwaysInfo{p}.GenesENSG,'stable');
        AllPathwaysInfo{p}.CoveredENSGgenes = a;
        AllPathwaysInfo{p}.ENSGSamplesCoverage = CommGeneFPKMsMat(:,c);
    end
    
    save('Results\Pathways\AllPathwaysInfo.mat','AllPathwaysInfo');
end

%% EXtra specific Regulations

RegulatedSamples(end+1).Name = 'GFP UpReg';
tmpindex = RegulatedSamples(2).Folds > 0;
RegulatedSamples(end).GeneIDs = RegulatedSamples(2).GeneIDs(tmpindex);
RegulatedSamples(end).FPKMsS1 = RegulatedSamples(2).FPKMsS1(tmpindex);
RegulatedSamples(end).FPKMsS2 = RegulatedSamples(2).FPKMsS2(tmpindex);
RegulatedSamples(end).Folds = RegulatedSamples(2).Folds(tmpindex);
RegulatedSamples(end).FPKMSFoldsMat = RegulatedSamples(2).FPKMSFoldsMat(tmpindex,:);

RegulatedSamples(end+1).Name = 'GFP DownReg';
tmpindex = RegulatedSamples(2).Folds < 0;
RegulatedSamples(end).GeneIDs = RegulatedSamples(2).GeneIDs(tmpindex);
RegulatedSamples(end).FPKMsS1 = RegulatedSamples(2).FPKMsS1(tmpindex);
RegulatedSamples(end).FPKMsS2 = RegulatedSamples(2).FPKMsS2(tmpindex);
RegulatedSamples(end).Folds = RegulatedSamples(2).Folds(tmpindex);
RegulatedSamples(end).FPKMSFoldsMat = RegulatedSamples(2).FPKMSFoldsMat(tmpindex,:);

RegulatedSamples(end+1).Name = 'RGpos UpReg';
tmpindex = RegulatedSamples(1).Folds > 0;
RegulatedSamples(end).GeneIDs = RegulatedSamples(1).GeneIDs(tmpindex);
RegulatedSamples(end).FPKMsS1 = RegulatedSamples(1).FPKMsS1(tmpindex);
RegulatedSamples(end).FPKMsS2 = RegulatedSamples(1).FPKMsS2(tmpindex);
RegulatedSamples(end).Folds = RegulatedSamples(1).Folds(tmpindex);
RegulatedSamples(end).FPKMSFoldsMat = RegulatedSamples(1).FPKMSFoldsMat(tmpindex,:);

RegulatedSamples(end+1).Name = 'RGpos DownReg';
tmpindex = RegulatedSamples(1).Folds < 0;
RegulatedSamples(end).GeneIDs = RegulatedSamples(1).GeneIDs(tmpindex);
RegulatedSamples(end).FPKMsS1 = RegulatedSamples(1).FPKMsS1(tmpindex);
RegulatedSamples(end).FPKMsS2 = RegulatedSamples(1).FPKMsS2(tmpindex);
RegulatedSamples(end).Folds = RegulatedSamples(1).Folds(tmpindex);
RegulatedSamples(end).FPKMSFoldsMat = RegulatedSamples(1).FPKMSFoldsMat(tmpindex,:);


for i = 1:length(RegulatedSamples)
    [~,tmpindx,~] = intersect(Sample(1).GenesENSG,RegulatedSamples(i).GeneIDs,'stable');
    RegulatedSamples(i).GeneShortIDs = Sample(1).GenesShortIDs(tmpindx);
end

%% Pathways outputs

[~,~,TmpPath] = xlsread('Pathways\Pathways-DAVID68.xlsx','GFP+');
TmpPath(1,end+1) = {'Entry'};
TmpPath(2:end,end) = extractBefore(TmpPath(2:end,2),':');
TmpPath(1,15:18) = {'Up-Reg ENSG','Down-Reg ENSG','Up-Reg ShortID','Down-Reg ShortID'};
TmpPath(1,19:22) = {'KEGG ID','Pathway Name','Class','Function'};
GeneGroups = TmpPath(2:end,6);
for i = 1:length(GeneGroups)
    tmp = strtrim(strsplit(GeneGroups{i},','));
    [~,~,c] = intersect(tmp,RegulatedSamples(3).GeneIDs,'stable');
    TmpPath(i+1,15) = {strjoin( RegulatedSamples(3).GeneIDs(c))};
    TmpPath(i+1,17) = {strjoin( RegulatedSamples(3).GeneShortIDs(c))};
    [~,~,c] = intersect(tmp,RegulatedSamples(4).GeneIDs,'stable');
    TmpPath(i+1,16) = {strjoin( RegulatedSamples(4).GeneIDs(c))};
    TmpPath(i+1,18) = {strjoin( RegulatedSamples(4).GeneShortIDs(c))};
    TmpPath(i+1,1) = {['= HYPERLINK("https://www.kegg.jp/dbget-bin/www_bget?pathway+' TmpPath{i+1,14} '","' TmpPath{i+1,14} '")']};
end
for i = 2:size(TmpPath,1)
    for p = 1:length(AllPathwaysInfo)     
        if strcmp(TmpPath{i,14},AllPathwaysInfo{p}.Entry) 
            TmpPath{i,19} = ['= HYPERLINK("https://www.kegg.jp/dbget-bin/www_bget?pathway+' AllPathwaysInfo{p}.Entry '","' AllPathwaysInfo{p}.Entry '"'];
            TmpPath{i,20} = AllPathwaysInfo{p}.Name;
            TmpPath{i,21} = AllPathwaysInfo{p}.MainClass;
            TmpPath{i,22} = AllPathwaysInfo{p}.Class;
        end
    end
end
GFP = TmpPath;

xlswrite('Pathways\Pathways-Res.xlsx',GFP,'GFP','A1');


[~,~,TmpPath] = xlsread('Pathways\Pathways-DAVID68.xlsx','RG+');
TmpPath(1,end+1) = {'Entry'};
TmpPath(2:end,end) = extractBefore(TmpPath(2:end,2),':');
TmpPath(1,15:18) = {'Up-Reg ENSG','Down-Reg ENSG','Up-Reg ShortID','Down-Reg ShortID'};
TmpPath(1,19:22) = {'KEGG ID','Pathway Name','Class','Function'};
GeneGroups = TmpPath(2:end,6);
for i = 1:length(GeneGroups)
    tmp = strtrim(strsplit(GeneGroups{i},','));
    [~,~,c] = intersect(tmp,RegulatedSamples(5).GeneIDs,'stable');
    TmpPath(i+1,15) = {strjoin( RegulatedSamples(5).GeneIDs(c))};
    TmpPath(i+1,17) = {strjoin( RegulatedSamples(5).GeneShortIDs(c))};
    [~,~,c] = intersect(tmp,RegulatedSamples(6).GeneIDs,'stable');
    TmpPath(i+1,16) = {strjoin( RegulatedSamples(6).GeneIDs(c))};
    TmpPath(i+1,18) = {strjoin( RegulatedSamples(6).GeneShortIDs(c))};
    TmpPath(i+1,1) = {['= HYPERLINK("https://www.kegg.jp/dbget-bin/www_bget?pathway+' TmpPath{i+1,14} '","' TmpPath{i+1,14} '")']};
end
for i = 2:size(TmpPath,1)
    for p = 1:length(AllPathwaysInfo)     
        if strcmp(TmpPath{i,14},AllPathwaysInfo{p}.Entry) 
            TmpPath{i,19} = ['= HYPERLINK("https://www.kegg.jp/dbget-bin/www_bget?pathway+' AllPathwaysInfo{p}.Entry '","' AllPathwaysInfo{p}.Entry '"'];
            TmpPath{i,20} = AllPathwaysInfo{p}.Name;
            TmpPath{i,21} = AllPathwaysInfo{p}.MainClass;
            TmpPath{i,22} = AllPathwaysInfo{p}.Class;
            
            %%----- input pathway genes => Output = foldchange matrices
            tmpsplited = strtrim(split(TmpPath{i,6},','));
            [~,indxtable,~] = intersect(table2array(T(:,1)),tmpsplited,'stable');
            xlswrite('Pathways\Pathways-Res.xlsx',table2cell(T(indxtable,:)),AllPathwaysInfo{p}.Entry,'D2');
        end
    end
    
    
end
RGpos = TmpPath;

xlswrite('Pathways\Pathways-Res.xlsx',RGpos,'RGpos','A1');



%%

% [a1,b,c] = intersect(OthersSigSamples(1).GeneIDs,RegulatedSamples(3).GeneIDs,'stable');
% [a2,b,c] = intersect(OthersSigSamples(3).GeneIDs,RegulatedSamples(3).GeneIDs,'stable');
% [a3,b,c] = intersect(OthersSigSamples(6).GeneIDs,RegulatedSamples(3).GeneIDs,'stable');
% [a4,b,c] = intersect(a1,a2,'stable');
% [a5,b,c] = intersect(a4,a3,'stable');
% 
% 
% [~,Other1,c] = intersect(OthersSigSamples(1).GeneIDs,a5,'stable');
% [~,Other2,c] = intersect(OthersSigSamples(3).GeneIDs,a5,'stable');
% [~,Other3,c] = intersect(OthersSigSamples(6).GeneIDs,a5,'stable');
% [aa,Ours1,c] = intersect(RegulatedSamples(1).GeneIDs,a5,'stable');
% [aa,Ours2,c] = intersect(RegulatedSamples(2).GeneIDs,a5,'stable');
% [~,Ours3,c] = intersect(RegulatedSamples(3).GeneIDs,a5,'stable');
% 
% 
% NEGwOthersFC = [];
% NEGwOthersFC = [[OthersSigSamples(6).Folds(Other3)]'];
% NEGwOthersFC = [NEGwOthersFC [OthersSigSamples(3).Folds(Other2)]'];
% NEGwOthersFC = [NEGwOthersFC [OthersSigSamples(1).Folds{Other1}]']; 
% NEGwOthersFC = [NEGwOthersFC [RegulatedSamples(1).Folds(Ours1)]];
% 
% GFPwOthersFC = [];
% GFPwOthersFC = [[OthersSigSamples(6).Folds(Other3)]'];
% GFPwOthersFC = [GFPwOthersFC [OthersSigSamples(3).Folds(Other2)]'];
% GFPwOthersFC = [GFPwOthersFC [OthersSigSamples(1).Folds{Other1}]']; 
% GFPwOthersFC = [GFPwOthersFC [RegulatedSamples(2).Folds(Ours2) 0]'];
% GFPwOthersFC(2,:) = [];
% 
% RGposwOthersFC = [];
% RGposwOthersFC = [[OthersSigSamples(6).Folds(Other3)]'];
% RGposwOthersFC = [RGposwOthersFC [OthersSigSamples(3).Folds(Other2)]'];
% RGposwOthersFC = [RGposwOthersFC [OthersSigSamples(1).Folds{Other1}]']; 
% RGposwOthersFC = [RGposwOthersFC [RegulatedSamples(3).Folds(Ours3)]];
% 
% 
% 
% 
% [a13,b,c] = intersect(OthersSigSamples(1).GeneIDs,OthersSigSamples(3).GeneIDs,'stable');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                    
%%
disp(' --------- Done !!! -------- ');



%% show the heatmaps of the Exp levels,



%% ------------------Functions-------------
function [Genes] = GroupingGenes (Uniqgenes, GeneFPKMsMat, FPKM_Thresh, NoiseRatio )

    FPKMsValuesMat = GeneFPKMsMat';

    %NC Genes
    NoiseIndex = sum(FPKMsValuesMat(:,1)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,1),2)*NoiseRatio);
    Genes.NC = Uniqgenes(not(NoiseIndex));

    %GFP Genes
    NoiseIndex = sum(FPKMsValuesMat(:,2)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,2),2)*NoiseRatio);
    Genes.GFP = Uniqgenes(not(NoiseIndex));

    %RGpos Genes
    NoiseIndex = sum(FPKMsValuesMat(:,3)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,3),2)*NoiseRatio);
    Genes.RGpos = Uniqgenes(not(NoiseIndex));

    Genes.GFP_RGpos = Genes.GFP(ismember(Genes.GFP,Genes.RGpos));
    Genes.NC_GFP = Genes.NC(ismember(Genes.NC,Genes.GFP));
    Genes.NC_RGpos = Genes.NC(ismember(Genes.NC,Genes.RGpos));

    Genes.RGpos_GFP_NC = Genes.RGpos((ismember(Genes.RGpos,Genes.GFP(ismember(Genes.GFP,Genes.NC)))));

    Genes.Just_NC = Genes.NC(~ismember(Genes.NC,Genes.GFP) & ~ismember(Genes.NC,Genes.RGpos));
    Genes.Just_GFP = Genes.GFP(~ismember(Genes.GFP,Genes.NC) & ~ismember(Genes.GFP,Genes.RGpos));
    Genes.Just_RGpos = Genes.RGpos(~ismember(Genes.RGpos,Genes.NC) & ~ismember(Genes.RGpos,Genes.GFP));

end


function [Genes] = GroupingGenes2 (Uniqgenes, GeneFPKMsMat, FPKM_Thresh, NoiseRatio )

    FPKMsValuesMat = GeneFPKMsMat';

    %NC Genes
    NoiseIndex = sum(FPKMsValuesMat(:,1)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,1),2)*NoiseRatio);
    Genes.NC = Uniqgenes(not(NoiseIndex));

    %NEG Genes
    NoiseIndex = sum(FPKMsValuesMat(:,2)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,2),2)*NoiseRatio);
    Genes.NEG = Uniqgenes(not(NoiseIndex));

    %GFP Genes
    NoiseIndex = sum(FPKMsValuesMat(:,3)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,3),2)*NoiseRatio);
    Genes.GFP = Uniqgenes(not(NoiseIndex));

    %RGpos Genes
    NoiseIndex = sum(FPKMsValuesMat(:,4)<FPKM_Thresh,2) >= round(size(FPKMsValuesMat(:,4),2)*NoiseRatio);
    Genes.RGpos = Uniqgenes(not(NoiseIndex));

    Genes.GFP_RGpos = Genes.GFP(ismember(Genes.GFP,Genes.RGpos));
    Genes.NEG_GFP = Genes.NEG(ismember(Genes.NEG,Genes.GFP));
    Genes.NEG_RGpos = Genes.NEG(ismember(Genes.NEG,Genes.RGpos));
    Genes.NC_NEG = Genes.NC(ismember(Genes.NC,Genes.NEG));
    Genes.NC_GFP = Genes.NC(ismember(Genes.NC,Genes.GFP));
    Genes.NC_RGpos = Genes.NC(ismember(Genes.NC,Genes.RGpos));

    Genes.RGpos_GFP_NC = Genes.RGpos((ismember(Genes.RGpos,Genes.GFP(ismember(Genes.GFP,Genes.NC)))));
    Genes.RGpos_GFP_NEG = Genes.RGpos((ismember(Genes.RGpos,Genes.GFP(ismember(Genes.GFP,Genes.NEG)))));
    Genes.RGpos_NEG_NC = Genes.RGpos((ismember(Genes.RGpos,Genes.NEG(ismember(Genes.NEG,Genes.NC)))));
    Genes.GFP_NEG_NC = Genes.GFP((ismember(Genes.GFP,Genes.NEG(ismember(Genes.NEG,Genes.NC)))));

    Genes.RGpos_GFP_NEG_NC = Genes.RGpos((ismember(Genes.RGpos,Genes.GFP(ismember(Genes.GFP,Genes.NC(ismember(Genes.NC,Genes.NEG)))))));

    Genes.Just_NC = Genes.NC(~ismember(Genes.NC,Genes.NEG) & ~ismember(Genes.NC,Genes.GFP) & ~ismember(Genes.NC,Genes.RGpos));
    Genes.Just_NEG = Genes.NEG(~ismember(Genes.NEG,Genes.NC) & ~ismember(Genes.NEG,Genes.GFP) & ~ismember(Genes.NEG,Genes.RGpos));
    Genes.Just_GFP = Genes.GFP(~ismember(Genes.GFP,Genes.NC) & ~ismember(Genes.GFP,Genes.NEG) & ~ismember(Genes.GFP,Genes.RGpos));
    Genes.Just_RGpos = Genes.RGpos(~ismember(Genes.RGpos,Genes.NC) & ~ismember(Genes.RGpos,Genes.NEG) & ~ismember(Genes.RGpos,Genes.GFP));

end

function [] = SVDs (GeneFPKMsMat, GroupName, Dimension, Path, Visibility)

    c=zeros(length(GroupName),3);
    c([1], :) = ones(1,3).*[0,0,1]; % NC
    c([2], :) = ones(1,3).*[0,0,0]; % NEG
    c([3], :) = ones(1,3).*[0,1,0]; % GFP
    c([4], :) = ones(1,3).*[1,0,0]; % RGpos

    [U,S,V] = svd(GeneFPKMsMat);
    
    switch Dimension
        case '2D'
            %--2D
            SVDrec = U(:,1:2)*S(1:2,1:2)*V(1:2,1:2)';
            figure('Visible',Visibility,'Position', get(0, 'Screensize'))
            % scatter(SVDrec(:,1),SVDrec(:,2),25)
            scatter(SVDrec(:,1),SVDrec(:,2),25,c,'filled')
            text(SVDrec(:,1),SVDrec(:,2),GroupName)
            title('SVD')
            xlabel({'Deviation Percentage : ';num2str((S(1,1)/sum(S(:)))*100)})
            ylabel({'Deviation Percentage : ';num2str((S(2,2)/sum(S(:)))*100)})
            saveas(gcf,[Path 'SVD- 2D .jpg']);
            saveas(gcf,[Path 'SVD- 2D .fig']);
            SavePDF(gcf,[Path 'SVD- 2D']);
            close all force;
        case '3D'
            %--3D
            SVDrec = U(:,1:3)*S(1:3,1:3)*V(1:3,1:3)';
            figure('Visible',Visibility,'Position', get(0, 'Screensize'))
            scatter3(SVDrec(:,1),SVDrec(:,2),SVDrec(:,3),25,c,'filled')
            text(SVDrec(:,1),SVDrec(:,2),SVDrec(:,3),GroupName)
            title('SVD')
            xlabel({'Deviation Percentage : ';num2str((S(1,1)/sum(S(:)))*100)})
            ylabel({'Deviation Percentage : ';num2str((S(2,2)/sum(S(:)))*100)})
            zlabel({'Deviation Percentage : ';num2str((S(3,3)/sum(S(:)))*100)})
            saveas(gcf,[Path 'SVD- 3D .jpg']);
            saveas(gcf,[Path 'SVD- 3D .fig']);
            SavePDF(gcf,[Path 'SVD- 3D']);
            close all force;
        otherwise
                disp('-----------------------ERROR---------------------------');
                error('--------Choose the correct dimensions (2D/3D)----------');
    end

end

function [] = BoxPlots (GeneFPKMsMat, GroupName, Path, Visibility)

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    boxplot(log2(GeneFPKMsMat+1)','Notch','on')
    set(gca,'XTickLabel',GroupName,'XTickLabelRotation',90);
    title('Gene Expressions Boxplot- All');
    xlabel('Samples')
    ylabel('log2(FPKMs)')
    saveas(gcf,[Path 'Box Plot- All.jpg']);
    saveas(gcf,[Path 'Box Plot- All.fig']);
    SavePDF(gcf,[Path 'Box Plot- All']);
    close all force;

end

function [] = ClusterGrams (GeneFPKMsMat, GroupName, GeneIDs, Type, Path, Titles)
    if Type == 1
        ClustObj = clustergram(GeneFPKMsMat','Standardize','Row','ColumnLabels',GroupName,'Linkage','complete','Dendrogram',0,'Colormap',redbluecmap);
    else
        ClustObj = clustergram(GeneFPKMsMat','Standardize','Row','ColumnLabels',GroupName,'RowLabels',GeneIDs,'Linkage','complete','Dendrogram',0,'Colormap',redbluecmap);
    end
    addTitle (ClustObj,{'All', Titles});
    addXLabel (ClustObj,'Samples');	
    addYLabel (ClustObj,'Genes');
    close all
    plot(ClustObj);
    set(gcf,'Position', get(0, 'Screensize'))
    % ax = ClustObj.plot;
    % colorbar('Peer', ax);
    saveas(gcf,[Path 'Clustergram-All' Titles '.jpg']);
    saveas(gcf,[Path 'Clustergram-All' Titles '.fig']);
    SavePDF(gcf,[Path 'Clustergram-All' Titles]);
    close all force;
end

function [] = GroupsGenesEXCEL (Genes, name, Titles ,Path, Visibility)
    GeneGroups = fieldnames(Genes);
    Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    filename = [Path name '.xlsx'];
    sheet = 1;
    xlRange = 'A1';
    xlswrite(filename,GeneGroups',sheet,xlRange)
    ExpGenesNum = [];
    for i = 1:length(GeneGroups)
        eval([' GenesName = Genes.' cell2mat(GeneGroups(i)) ';']);
        if i<=26
            xlRange = [Alphabet(i) int2str(2)];
        else
            xlRange = [Alphabet(floor(i/26)) Alphabet(i-floor(i/26)*26+(floor(i/26)-1)) int2str(2)];
        end
        xlswrite(filename,length(GenesName),sheet,xlRange)
        ExpGenesNum = [ExpGenesNum length(GenesName)];
        if isempty(GenesName)
            continue;
        end
        if i<=26
            xlRange = [Alphabet(i) int2str(3)];
        else
            xlRange = [Alphabet(floor(i/26)) Alphabet(i-floor(i/26)*26+(floor(i/26)-1)) int2str(3)];
        end
        xlswrite(filename,GenesName,sheet,xlRange)
    end
    
    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    bar(ExpGenesNum,0.4)
    title({Titles;'All Samples'});
    set(gca,'XTick',[1:length(ExpGenesNum)],'XTickLabel',GeneGroups,'XTickLabelRotation',90);
    xlabel('Samples')
    ylabel('# Expressed gene')
    saveas(gcf,[Path 'Num Expressed- ' Titles '.jpg']);
    saveas(gcf,[Path 'Num Expressed- ' Titles '.fig']);
    SavePDF(gcf,[Path 'Num Expressed- ' Titles ]);
    close all force;

end

function [] = GroupSpecExpGenes(tmpGeneFPKMs, GeneENSG, tmpL, Type, Path, Visibility)

    [tmpGeneFPKMs,Ind] = sort(tmpGeneFPKMs,'descend') ;

    if ~isempty(Ind)
        GeneENSG = GeneENSG(Ind);
        
        figure('Visible',Visibility,'Position', get(0, 'Screensize'))
        bar(tmpGeneFPKMs)
        title([tmpL ' - EXP Genes']);
        if Type ~= 1
            set(gca,'XTick',[1:length(GeneENSG)],'XTickLabel',GeneENSG,'XTickLabelRotation',90);
        end
        xlabel('Genes')
        ylabel('FPKMs')
        saveas(gcf,[Path tmpL ' - EXP Genes.jpg']);
        saveas(gcf,[Path tmpL ' - EXP Genes.fig']);
        SavePDF(gcf,[Path tmpL ' - EXP Genes']);
        close all force;

        N = min(50,length(tmpGeneFPKMs));
        figure('Visible',Visibility,'Position', get(0, 'Screensize'))
        bar(tmpGeneFPKMs(1:N))
        title({[tmpL ' - EXP Genes'],'Top 50'});
        set(gca,'XTick',[1:length(GeneENSG(1:N))],'XTickLabel',GeneENSG(1:N),'XTickLabelRotation',90);
        xlabel('Genes')
        ylabel('FPKMs')
        saveas(gcf,[Path tmpL ' - EXP Genes - Top Exp.jpg']);
        saveas(gcf,[Path tmpL ' - EXP Genes - Top Exp.fig']);
        SavePDF(gcf,[Path tmpL ' - EXP Genes - Top Exp']);
        close all force;
    end


end

function [RegulatedFoldsValuesMat] = RegulatedFolds (RegulatedGenes, RegulatedSamples, Samps )
    RegulatedFoldsValuesMat = zeros(length(RegulatedGenes),length(Samps));
    for i = 1:length(Samps)
        tmpindex = find(ismember(RegulatedSamples(Samps(i)).GeneIDs,RegulatedGenes));
        RegulatedFoldsValuesMat(:,i) = RegulatedSamples(Samps(i)).Folds(tmpindex);
    end
    RegulatedFoldsValuesMat(RegulatedFoldsValuesMat > 20) = 20;
    RegulatedFoldsValuesMat(RegulatedFoldsValuesMat < -20) = -20;
end

function [RegulatedGenes] = RegulationBars(RegulatedGenes, RegulatedFoldsValuesMat, CommonRegulatedGenes, Colors, Titles, Legends, Figurename, Path, Visibility)

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(RegulatedFoldsValuesMat);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(RegulatedFoldsValuesMat)],'XTickLabel',CommonRegulatedGenes,'XTickLabelRotation',90);
    title(Titles);
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes Regulations.jpg']);
    saveas(gcf,[Path Figurename 'Genes Regulations.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes Regulations']);
    close all force;
    
      
    %----- UP REGULATION
    tmpLogFoldsUp = sign(RegulatedFoldsValuesMat);
    tmpIndexUp = sum(tmpLogFoldsUp,2)==length(Colors);
    tmpLogFoldsUp = RegulatedFoldsValuesMat(tmpIndexUp,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsUp);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsUp)],'XTickLabel',CommonRegulatedGenes(tmpIndexUp),'XTickLabelRotation',90);
    title({Titles; 'Up Regulation'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes UP Regulation.jpg']);
    saveas(gcf,[Path Figurename 'Genes UP Regulation.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes UP Regulation']);
    close all force;
    
    eval([' RegulatedGenes.' Titles '_UpReg = CommonRegulatedGenes(tmpIndexUp);']);
    
    
    tmpIndexUp_NoNC = sum(RegulatedFoldsValuesMat==20,2)==length(Colors);
    tmpLogFoldsUp_NoNC = RegulatedFoldsValuesMat(tmpIndexUp_NoNC,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsUp_NoNC);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsUp_NoNC)],'XTickLabel',CommonRegulatedGenes(tmpIndexUp_NoNC),'XTickLabelRotation',90);
    title({Titles; 'Genes UP Regulation';'No Expression in NC'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes UP Regulation - No NC.jpg']);
    saveas(gcf,[Path Figurename 'Genes UP Regulation - No NC.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes UP Regulation - No NC']);
    close all force;
    
    eval([' RegulatedGenes.' Titles '_UpReg_NoNC = CommonRegulatedGenes(tmpIndexUp_NoNC);']);
    
    tmpIndexUp_WithNC = sum(RegulatedFoldsValuesMat>0 & RegulatedFoldsValuesMat<20,2)>0;
    tmpLogFoldsUp_WithNC = RegulatedFoldsValuesMat(tmpIndexUp_WithNC,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsUp_WithNC);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsUp_WithNC)],'XTickLabel',CommonRegulatedGenes(tmpIndexUp_WithNC),'XTickLabelRotation',90);
    title({Titles; 'Genes UP Regulation';'With Expression in NC'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes UP Regulation - With NC.jpg']);
    saveas(gcf,[Path Figurename 'Genes UP Regulation - With NC.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes UP Regulation - With NC']);
    close all force;

    eval([' RegulatedGenes.' Titles '_UpReg_WithNC = CommonRegulatedGenes(tmpIndexUp_WithNC);']);
    
    %---Down Regulation

    tmpLogFoldsDown = sign(RegulatedFoldsValuesMat);
    tmpIndexDown = sum(tmpLogFoldsDown,2)==-length(Colors);
    tmpLogFoldsDown = RegulatedFoldsValuesMat(tmpIndexDown,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsDown);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsDown)],'XTickLabel',CommonRegulatedGenes(tmpIndexDown),'XTickLabelRotation',90);
    title({Titles; 'Genes Down Regulation'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes Down Regulation.jpg']);
    saveas(gcf,[Path Figurename 'Genes Down Regulation.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes Down Regulation']);
    close all force;

    eval([' RegulatedGenes.' Titles '_DownReg = CommonRegulatedGenes(tmpIndexDown);']);
    
    tmpIndexDown_JustNC = sum(RegulatedFoldsValuesMat==-20,2)==length(Colors);
    tmpLogFoldsDown_JustNC = RegulatedFoldsValuesMat(tmpIndexDown_JustNC,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsDown_JustNC);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsDown_JustNC)],'XTickLabel',CommonRegulatedGenes(tmpIndexDown_JustNC),'XTickLabelRotation',90);
    title({Titles; 'Genes Down Regulation';'Just Expressed in NC'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes Down Regulation - Just NC.jpg']);
    saveas(gcf,[Path Figurename 'Genes Down Regulation - Just NC.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes Down Regulation - Just NC']);
    close all force;

    eval([' RegulatedGenes.' Titles '_DownReg_NoExp = CommonRegulatedGenes(tmpIndexDown_JustNC);']);

    tmpIndexDown_WithOthers = sum(RegulatedFoldsValuesMat<0 & RegulatedFoldsValuesMat>-20,2)>0;
    tmpLogFoldsDown_WithOthers = RegulatedFoldsValuesMat(tmpIndexDown_WithOthers,:);

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    hb = bar(tmpLogFoldsDown_WithOthers);
    for c = 1:length(Colors)
        set(hb(c), 'FaceColor',Colors(c))
    end
    set(gca,'XTick',[1:length(tmpLogFoldsDown_WithOthers)],'XTickLabel',CommonRegulatedGenes(tmpIndexDown_WithOthers),'XTickLabelRotation',90);
    title({Titles;'Genes Down Regulation';'At least Expressed in one sample'});
    legend(Legends)
    xlabel('Genes')
    ylabel('log2 genes regulation')
    saveas(gcf,[Path Figurename 'Genes Down Regulation - non-zero Expresion.jpg']);
    saveas(gcf,[Path Figurename 'Genes Down Regulation - non-zero Expresion.jpg']);
    SavePDF(gcf,[Path Figurename 'Genes Down Regulation - non-zero Expresion']);
    close all force;
    
    eval([' RegulatedGenes.' Titles '_DownReg_WithExp = CommonRegulatedGenes(tmpIndexDown_WithOthers);']);


end

function [] = TopNexpressed(Sample, SampInd, N, XlsxFileName, tmpL, Path, Visibility)
    
    N = min(N,length(Sample(SampInd(1)).GenesFPKMs));
    TopsIndx = zeros(N,length(SampInd));
    TopsGenes = {};
    for i = 1:length(SampInd)
        [~,Ind] = sort(Sample(SampInd(i)).GenesFPKMs,'descend');
        TopsIndx(:,SampInd(i)) = Ind(1:N);
%         TopsGenes = union(TopsGenes,Sample(SampInd(i)).GenesIDs);
        TopsGenes = union(TopsGenes,Sample(SampInd(i)).GenesIDs(Ind(1:N)));
    end
    
    TopFPKMsMat = zeros(length(TopsGenes),length(SampInd));
    
    for i = 1:length(SampInd)
        [~,~,RefB] = intersect(TopsGenes,Sample(SampInd(i)).GenesIDs,'stable');
        TopFPKMsMat(:,i) = Sample(SampInd(i)).GenesFPKMs(RefB);
    end
    
    [~,~,RefB] = intersect(TopsGenes,Sample(SampInd(i)).GenesIDs,'stable');
    
    TopsGenes = Sample(SampInd(1)).GenesShortIDs(RefB);
 
    filename = [Path XlsxFileName];
    sheet = 1;
    xlswrite(filename,TopsGenes,sheet,'A2')
    xlswrite(filename,{'Genes ID' , Sample(SampInd).Name},sheet,'A1')
    xlswrite(filename,TopFPKMsMat,sheet,'B2')
        

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    bar(TopFPKMsMat)
    title([tmpL ' - EXP Genes']);
    set(gca,'XTick',[1:length(TopsGenes)],'XTickLabel',TopsGenes,'XTickLabelRotation',90);
    xlabel('Genes')
    ylabel('FPKMs')
    legend({Sample(SampInd).Name});
    saveas(gcf,[Path tmpL ' - EXP Genes.jpg']);
    saveas(gcf,[Path tmpL ' - EXP Genes.fig']);
    SavePDF(gcf,[Path tmpL ' - EXP Genes']);
    close all force;

    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    bar(log2(TopFPKMsMat+1))
    title([tmpL ' - EXP Genes']);
    set(gca,'XTick',[1:length(TopsGenes)],'XTickLabel',TopsGenes,'XTickLabelRotation',90);
    xlabel('Genes')
    ylabel('log2(FPKMs)')
    legend({Sample(SampInd).Name});
    saveas(gcf,[Path tmpL ' - EXP Genes - log2.jpg']);
    saveas(gcf,[Path tmpL ' - EXP Genes - log2.fig']);
    SavePDF(gcf,[Path tmpL ' - EXP Genes - log2']);
    close all force;
    
    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    h = heatmap(TopsGenes,{Sample(SampInd).Name},TopFPKMsMat');
    h.ColorScaling = 'scaledcolumns';
    h.Title = ['Top Expressed genes'];
    h.YLabel = 'Samples';
    h.XLabel = 'Genes';
    h.Colormap = parula;
    colorbar;
    saveas(gcf,[Path tmpL 'Heatmap all fpkms.jpg']);
    saveas(gcf,[Path tmpL 'Heatmap all fpkms.fig']);
    SavePDF(gcf,[Path tmpL 'Heatmap all fpkms']);
    close all force;
    
    figure('Visible',Visibility,'Position', get(0, 'Screensize'))
    h = heatmap(TopsGenes,{Sample(SampInd).Name},log2(TopFPKMsMat+1)');
    h.Title = ['Log2 of Top Expressed genes'];
    h.YLabel = 'Samples';
    h.XLabel = 'Genes';
    h.Colormap = parula;
    saveas(gcf,[Path tmpL 'Heatmap all fpkms- log2.jpg']);
    saveas(gcf,[Path tmpL 'Heatmap all fpkms- log2.fig']);
    SavePDF(gcf,[Path tmpL 'Heatmap all fpkms- log2']);
    close all force;
end

function [] = SavePDF(h, Filename)
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(h,Filename,'-dpdf','-r0')
end

function [] = AllGenesExpressionEXCEL (Sample, FPKMsMat, name, Path)
    Header = {'Genes ID'};
    for i = 1:length(Sample)
        Header{end+1} = Sample(i).Name;
    end
    
    filename = [Path name '.xlsx'];
    
    sheet1 = 'ENSG';
    sheet2 = 'Gene Symbol';
    
    xlRange = 'A1';
    xlswrite(filename,Header,sheet1,xlRange);
    xlswrite(filename,Header,sheet2,xlRange);
    
    xlRange = 'A2';
    xlswrite(filename,Sample(1).GenesENSG,sheet1,xlRange);
    xlswrite(filename,Sample(1).GenesShortIDs,sheet2,xlRange);
    
    xlRange = 'B2';
    xlswrite(filename,FPKMsMat',sheet1,xlRange);
    xlswrite(filename,FPKMsMat',sheet2,xlRange);
end


%% 
