clc
clear all
clc
biopsy_patches_directory = 'E:\Cache_Tiles'; % path to biopsy patches directories
dirim = dir(biopsy_patches_directory); % Get directory information (all biopsy patch directories)

filename = 'Cache_AllFeatures.xlsx';


% Initializng the structure with zeros
img = imread('D:\25A_8.png'); % An image with zero mask
mask=imread('D:\25A_8_mask.png');
    
bounds=DLMask2bounds_function_noorginput(mask);

% lymphocyte connection visulization
r=0.10;
alpha=0.43;
CC = bwconncomp(mask);
[VY,VX,x1,y1,edges] = construct_ccgs_optimized(bounds,alpha, r);

edges= edges | edges';

[numcomp,group] = graphconncomp(sparse(edges));%MyoFociGraphCount


Name='';
P_Name='';
id='';

FociCount=0; 
MyoFociCount=0;
FociGraphCount=0;
Myo_FociGraphCount=0;
RawMyo_FociGraphCount=0;
Tissue_MyoFociGraphCount=0;
LymphFreeTissue_FociGraphCount=0;
Title_MyoFociGraphCount=0;
Lymph_MyoFociGraphCount=0;

LymphAreaRatio1=0;
MyoAreaRatio1=0;
MyoAreaRatio2=0;
MyoAreaRatio3=0;
RawMyoAreaRatio1=0;
RawMyoAreaRatio2=0;
RawMyoAreaRatio3=0;
FociArea=0;
FociAreaRatio1=0;
FociAreaRatio2=0;
FociAreaRatio3=0;
FociAreaRatio4=0;
FociAreaRatio5=0;
FociAreaRatio6=0;
MyoFociArea=0;
MyoFociAreaRatio1=0;
MyoFociAreaRatio2=0;
MyoFociAreaRatio3=0;
MyoFociAreaRatio4=0;
MyoFociAreaRatio5=0;
MyoArea=0;
RawMyoArea=0;
LymphArea=0;
TissueArea=0;
LymphocyteArea=0;
LymphAreaRatio2=0;
LymphAreaRatio3=0;
LymphAreaRatio4=0;
LymphAreaRatio5=0;
LymphAreaRatio6=0;
LymphAreaRatio7=0;
LymphAreaRatio8=0;
LymphAreaRatio9=0;
LymphAreaRatio10=0;
LymphFreeArea=0;
TotalArea=0;



S= struct('patient_name',P_Name,'patch_name',Name,'patch_id',id,'MyoFociGraphCount',numcomp,'FociCount',FociCount,'MyoFociCount',MyoFociCount,'FociGraphCount',FociGraphCount,...
'Myo_FociGraphCount',Myo_FociGraphCount,'RawMyo_FociGraphCount',RawMyo_FociGraphCount,'Tissue_MyoFociGraphCount',Tissue_MyoFociGraphCount,'LymphFreeTissue_FociGraphCount',...
LymphFreeTissue_FociGraphCount,'Title_MyoFociGraphCount',Title_MyoFociGraphCount,'Lymph_MyoFociGraphCount',Lymph_MyoFociGraphCount,'LymphAreaRatio1',LymphAreaRatio1,'MyoAreaRatio1',...
MyoAreaRatio1,'MyoAreaRatio2',MyoAreaRatio2','MyoAreaRatio3',MyoAreaRatio3,'RawMyoAreaRatio1',RawMyoAreaRatio1,'RawMyoAreaRatio2',RawMyoAreaRatio2,...
'RawMyoAreaRatio3',RawMyoAreaRatio3,'FociArea',FociArea,'FociAreaRatio1',FociAreaRatio1,'FociAreaRatio2',FociAreaRatio2,'FociAreaRatio3',FociAreaRatio3,'FociAreaRatio4',FociAreaRatio4,...
'FociAreaRatio5',FociAreaRatio5,'FociAreaRatio6',FociAreaRatio6,'MyoFociArea',MyoFociArea,'MyoFociAreaRatio1',MyoFociAreaRatio1,...
'MyoFociAreaRatio2',MyoFociAreaRatio2,'MyoFociAreaRatio3',MyoFociAreaRatio3,'MyoFociAreaRatio4',MyoFociAreaRatio4,'MyoFociAreaRatio5',MyoFociAreaRatio5,'MyoArea',MyoArea,...
'RawMyoArea',RawMyoArea,'LymphArea',LymphArea,'TissueArea',TissueArea,'LymphocyteArea',LymphocyteArea,'LymphAreaRatio2',LymphAreaRatio2,'LymphAreaRatio3',LymphAreaRatio3,...
'LymphAreaRatio4',LymphAreaRatio4,'LymphAreaRatio5',LymphAreaRatio5,'LymphAreaRatio6',LymphAreaRatio6,'LymphAreaRatio7',LymphAreaRatio7,'LymphAreaRatio8',LymphAreaRatio8,...
'LymphAreaRatio9',LymphAreaRatio9,'LymphAreaRatio10',LymphAreaRatio10,'LymphFreeArea',LymphFreeArea,'TotalArea',TotalArea);




b=2;
for j = 1:length(dirim)
    fprintf('Creating mask of lymphocytes nuclei for %s\n',dirim(j).name);
    patient_folder = fullfile(biopsy_patches_directory,dirim(j).name); % actually the biopsy folder
    patches = dir([patient_folder '\*.png']); % Get all '.png' extensions
    for i = 1:length(patches)
        if length(patches(i).name)<3
            continue
        end
        if regexp(patches(i).name,'lymph_patch')
            continue
        end
        if regexp(patches(i).name,'_mask')
            continue
        end
        if regexp(patches(i).name,'_maask')
            continue
        end
        if regexp(patches(i).name,'MYOmask')
            continue
        end
        if regexp(patches(i).name,'RawMyoMask')
            continue
        end
        if regexp(patches(i).name,'_mask2')% This is actually redundant because we already took care of "_mask"
            continue
        end
        img_path = fullfile(patient_folder, patches(i).name);
        [~,patch_id,~] = fileparts(img_path);
        mask1_path = fullfile(patient_folder, [patch_id '_mask.png']);% On July 5th I changed "_mask.png" to "_mask2.png" to use 5K condition this round and neglect small foci
        mask2_path = fullfile(patient_folder, [patch_id '_MYOmask.png']);
        mask3_path = fullfile(patient_folder, [patch_id '_RawMyoMask.png']);
        mask4_path = fullfile(patient_folder, ['lymph_patch_3K_' patch_id '.png']);
        mask5_path = fullfile(patient_folder, [patch_id '_maask.png']);

        
        
        img = imread(img_path);
        mask1=imbinarize(imread(mask1_path));% Mask of lymphocyte foci
        mask2=imread(mask2_path);% mask of myocardium tissue
        mask=mask1.*mask2;% The mask for the foci in myocardium
        bounds=DLMask2bounds_function_noorginput(mask);
        
        testFlag=bwconncomp(mask);% Here we check for 2 things: 1)the foci "mask" is not empty so that we don't spend making the graph 2)if the foci "mask" is empty we don't care for myocardium "MYOmask" that is multiplied in it
        
        
        mask3=imread(mask3_path);
        mask4=imread(mask4_path);
        
        MyoArea=sum(sum(mask2)); %Area of myocardium in MYOmask
        RawMyoArea=sum(sum(mask3)); %Area of myocardium in RawMyoMask
        
        LymphMask = mask4(:,:,1);%Extracting red channel from the mask
        TissueMask = (mask4(:,:,3)<235)&(mask4(:,:,2)<230);% Extracting dark blue where is covered by the tissue
        LymphArea=sum(sum(imbinarize(LymphMask)));
        TissueArea=sum(sum(TissueMask));
        
        Lymph_free_TissueMask=TissueMask-(TissueMask.*imbinarize(LymphMask)); %Requested by Eliot
        LymphFreeArea=sum(sum(Lymph_free_TissueMask));
        
        LymphAreaRatio1=LymphArea/TissueArea; %Old version of lymphocyte area ratio(Thomas's)
        
        [aa bb cc]=size(img); %Getting the dimentions of the patch to calculate the whole area
        TotalArea=aa*bb;
        
        MyoAreaRatio1=MyoArea/TotalArea;
        MyoAreaRatio2=MyoArea/TissueArea;
        MyoAreaRatio3=MyoArea/LymphFreeArea;

        
        RawMyoAreaRatio1=RawMyoArea/TotalArea;
        RawMyoAreaRatio2=RawMyoArea/TissueArea;
        RawMyoAreaRatio3=RawMyoArea/LymphFreeArea;
        
        FociArea=sum(sum(mask1));
        FociAreaRatio1=FociArea/TotalArea;
        FociAreaRatio2=FociArea/TissueArea;
        FociAreaRatio3=FociArea/MyoArea;
        FociAreaRatio4=FociArea/RawMyoArea;
        FociAreaRatio5=FociArea/TotalArea;
        FociAreaRatio6=FociArea/LymphFreeArea;

        
        MyoFociArea=sum(sum(imbinarize(mask)));
        MyoFociAreaRatio1=MyoFociArea/MyoArea;
        MyoFociAreaRatio2=MyoFociArea/TissueArea;
        MyoFociAreaRatio3=MyoFociArea/RawMyoArea;
        MyoFociAreaRatio4=MyoFociArea/TotalArea;
        MyoFociAreaRatio5=MyoFociArea/LymphFreeArea;
        
        mask5=imread(mask5_path);

        LymphocyteArea=sum(sum(mask5));%Eliot mentioned as well
        LymphAreaRatio2=LymphocyteArea/TissueArea;%Eliot mentioned as well
        LymphAreaRatio3=LymphocyteArea/TotalArea;
        LymphAreaRatio4=LymphocyteArea/MyoArea;%Eliot mentioned as well
        LymphAreaRatio5=LymphocyteArea/RawMyoArea;%Eliot mentioned as well
        LymphAreaRatio6=LymphocyteArea/LymphFreeArea;% Eliot proposed
        
        %Based on LymphArea in blue masks
        %Old Version Calculated Above (LymphAreaRatio1)=LymphArea/TissueArea;
        LymphAreaRatio7=LymphArea/TotalArea;
        LymphAreaRatio8=LymphArea/MyoArea;
        LymphAreaRatio9=LymphArea/RawMyoArea;
        LymphAreaRatio10=LymphArea/LymphFreeArea;
        
        %mask4_temp=bwconncomp(mask4);%THIS IS WERE I EDITED on June 20th BASED ON ELIOT's COMMENT
        mask4_temp=bwconncomp(LymphMask);%MUST BE THE CORRECT ONE (It gettes the red regions out) 
        FociCount=mask4_temp.NumObjects;% number of components in old 3k masks
        
        se1 = strel('disk',10); % a structuring element for morphological dilation
        Dil_mask=imdilate(mask,se1);
        Dil_mask_temp=bwconncomp(Dil_mask);
        MyoFociCount=Dil_mask_temp.NumObjects;
        
        %Finding the number of proximal graphs in the "_mask" before
        %multiplying in MYOmask

        bounds=bwnuclei2bounds(logical(mask1));
        bounds=bounds';


        % lymphocyte connection visulization

        [VY,VX,x1,y1,edges] = construct_ccgs_optimized(bounds,alpha, r);
        edges= edges | edges';
        [FociGraphCount,group] = graphconncomp(sparse(edges));%FociGraphCount

             
        if testFlag.NumObjects~=0
            
            %Finding the number of proximal graphs in the "mult"
            bounds=bwnuclei2bounds(logical(mask));
            bounds=bounds';


            % lymphocyte connection visulization
            r=0.10;
            alpha=0.43;
            %CC = bwconncomp(mask);
            [VY,VX,x1,y1,edges] = construct_ccgs_optimized(bounds,alpha, r);
            edges= edges | edges';
            [numcomp,group] = graphconncomp(sparse(edges));%MyoFociGraphCount

            


            Myo_FociGraphCount=FociGraphCount/MyoArea;% Eliot proposed
            RawMyo_FociGraphCount=FociGraphCount/RawMyoArea;% Eliot proposed
            Tissue_FociGraphCount=FociGraphCount/TissueArea;% Eliot proposed
            LymphFreeTissue_FociGraphCount=FociGraphCount/LymphFreeArea;% Eliot proposed
            Tile_FociGraphCount=FociGraphCount/TotalArea;
            Lymph_FociGraphCount=FociGraphCount/LymphArea;
            
            
            
            Myo_MyoFociGraphCount=numcomp/MyoArea;
            RawMyo_MyoFociGraphCount=numcomp/RawMyoArea;
            Tissue_MyoFociGraphCount=numcomp/TissueArea;
            LymphFreeTissue_MyoFociGraphCount=numcomp/LymphFreeArea;
            Title_MyoFociGraphCount=numcomp/TotalArea;
            Lymph_MyoFociGraphCount=numcomp/LymphArea;
            
            
            P_Name=dirim(j).name;
            Name=patches(i).name;
            id=patch_id;
            
            
            S(1,b).patient_name=P_Name;
            S(1,b).patch_name=Name;
            S(1,b).patch_id=id;
            S(1,b).MyoFociGraphCount=numcomp;
            S(1,b).FociGraphCount=FociGraphCount;
            
            S(1,b).FociCount=FociCount;
            S(1,b).MyoFociCount=MyoFociCount;
            S(1,b).Myo_FociGraphCount=Myo_FociGraphCount;
            S(1,b).RawMyo_FociGraphCount=RawMyo_FociGraphCount;
            S(1,b).Tissue_MyoFociGraphCount=Tissue_MyoFociGraphCount;
            S(1,b).LymphFreeTissue_FociGraphCount=LymphFreeTissue_FociGraphCount;
            S(1,b).Title_MyoFociGraphCount=Title_MyoFociGraphCount;
            S(1,b).Lymph_MyoFociGraphCount=Lymph_MyoFociGraphCount;
            S(1,b).LymphAreaRatio1=LymphAreaRatio1;
            S(1,b).MyoAreaRatio1=MyoAreaRatio1;
            S(1,b).MyoAreaRatio2=MyoAreaRatio2;
            S(1,b).MyoAreaRatio3=MyoAreaRatio3;
            S(1,b).RawMyoAreaRatio1=RawMyoAreaRatio1;
            S(1,b).RawMyoAreaRatio2=RawMyoAreaRatio2;
            S(1,b).RawMyoAreaRatio3=RawMyoAreaRatio3;
            S(1,b).FociArea=FociArea;
            S(1,b).FociAreaRatio1=FociAreaRatio1;
            S(1,b).FociAreaRatio2=FociAreaRatio2;
            S(1,b).FociAreaRatio3=FociAreaRatio3;
            S(1,b).FociAreaRatio4=FociAreaRatio4;
            S(1,b).FociAreaRatio5=FociAreaRatio5;
            S(1,b).FociAreaRatio6=FociAreaRatio6;
   
            S(1,b).MyoFociArea=MyoFociArea;
            S(1,b).MyoFociAreaRatio1=MyoFociAreaRatio1;
            S(1,b).MyoFociAreaRatio2=MyoFociAreaRatio2;
            S(1,b).MyoFociAreaRatio3=MyoFociAreaRatio3;
            S(1,b).MyoFociAreaRatio4=MyoFociAreaRatio4;
            S(1,b).MyoFociAreaRatio5=MyoFociAreaRatio5;
            S(1,b).MyoArea=MyoArea;
            S(1,b).RawMyoArea=RawMyoArea;
            S(1,b).LymphArea=LymphArea;
            S(1,b).TissueArea=TissueArea;
            S(1,b).LymphocyteArea=LymphocyteArea;
            S(1,b).LymphAreaRatio2=LymphAreaRatio2;
            S(1,b).LymphAreaRatio3=LymphAreaRatio3;
            S(1,b).LymphAreaRatio4=LymphAreaRatio4;
            S(1,b).LymphAreaRatio5=LymphAreaRatio5;
            S(1,b).LymphAreaRatio6=LymphAreaRatio6;
            S(1,b).LymphAreaRatio7=LymphAreaRatio7;
            S(1,b).LymphAreaRatio8=LymphAreaRatio8;
            S(1,b).LymphAreaRatio9=LymphAreaRatio9;
            S(1,b).LymphAreaRatio10=LymphAreaRatio10;
            S(1,b).LymphFreeArea=LymphFreeArea;
            S(1,b).TotalArea=TotalArea;


            b=b+1;
        else 
            P_Name=dirim(j).name;
            Name=patches(i).name;
            id=patch_id;

            S(1,b).patient_name=P_Name;
            S(1,b).patch_name=Name;
            S(1,b).patch_id=id;
            S(1,b).MyoFociGraphCount=0;%numcomp
            S(1,b).FociGraphCount=0;%FociGraphCount

            S(1,b).FociCount=FociCount;
            S(1,b).MyoFociCount=MyoFociCount;
            S(1,b).Myo_FociGraphCount=Myo_FociGraphCount;
            S(1,b).RawMyo_FociGraphCount=RawMyo_FociGraphCount;
            S(1,b).Tissue_MyoFociGraphCount=Tissue_MyoFociGraphCount;
            S(1,b).LymphFreeTissue_FociGraphCount=LymphFreeTissue_FociGraphCount;
            S(1,b).Title_MyoFociGraphCount=Title_MyoFociGraphCount;
            S(1,b).Lymph_MyoFociGraphCount=Lymph_MyoFociGraphCount;
            S(1,b).LymphAreaRatio1=LymphAreaRatio1;
            S(1,b).MyoAreaRatio1=MyoAreaRatio1;
            S(1,b).MyoAreaRatio2=MyoAreaRatio2;
            S(1,b).MyoAreaRatio3=MyoAreaRatio3;
            S(1,b).RawMyoAreaRatio1=RawMyoAreaRatio1;
            S(1,b).RawMyoAreaRatio2=RawMyoAreaRatio2;
            S(1,b).RawMyoAreaRatio3=RawMyoAreaRatio3;
            S(1,b).FociArea=FociArea;
            S(1,b).FociAreaRatio1=FociAreaRatio1;
            S(1,b).FociAreaRatio2=FociAreaRatio2;
            S(1,b).FociAreaRatio3=FociAreaRatio3;
            S(1,b).FociAreaRatio4=FociAreaRatio4;
            S(1,b).FociAreaRatio5=FociAreaRatio5;
            S(1,b).FociAreaRatio6=FociAreaRatio6;
   
            S(1,b).MyoFociArea=MyoFociArea;
            S(1,b).MyoFociAreaRatio1=MyoFociAreaRatio1;
            S(1,b).MyoFociAreaRatio2=MyoFociAreaRatio2;
            S(1,b).MyoFociAreaRatio3=MyoFociAreaRatio3;
            S(1,b).MyoFociAreaRatio4=MyoFociAreaRatio4;
            S(1,b).MyoFociAreaRatio5=MyoFociAreaRatio5;
            S(1,b).MyoArea=MyoArea;
            S(1,b).RawMyoArea=RawMyoArea;
            S(1,b).LymphArea=LymphArea;
            S(1,b).TissueArea=TissueArea;
            S(1,b).LymphocyteArea=LymphocyteArea;
            S(1,b).LymphAreaRatio2=LymphAreaRatio2;
            S(1,b).LymphAreaRatio3=LymphAreaRatio3;
            S(1,b).LymphAreaRatio4=LymphAreaRatio4;
            S(1,b).LymphAreaRatio5=LymphAreaRatio5;
            S(1,b).LymphAreaRatio6=LymphAreaRatio6;
            S(1,b).LymphAreaRatio7=LymphAreaRatio7;
            S(1,b).LymphAreaRatio8=LymphAreaRatio8;
            S(1,b).LymphAreaRatio9=LymphAreaRatio9;
            S(1,b).LymphAreaRatio10=LymphAreaRatio10;
            S(1,b).LymphFreeArea=LymphFreeArea;
            S(1,b).TotalArea=TotalArea;            
            
            
            
            b=b+1;
        end
        
    end
end

T = struct2table(S);
save('Cache_AllFeatures.mat', 'T');
writetable(T,filename);

