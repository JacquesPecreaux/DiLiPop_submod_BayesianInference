% script to generate a template mat file from a template xmlfile

[xmlfile, p] = uigetfile('*.xml','Please choose a job file to process');
[~,xmlfile_bkp,~] = fileparts(xmlfile);
filexml = fullfile(p, xmlfile);

%% Normally, nothing to change below

docNode= xmlread(filexml);

%% get general_param, common to the project

work_path_List=docNode.getElementsByTagName('General_Params');
work_path_Node=work_path_List.item(work_path_List.getLength-1);
childNode = work_path_Node.getFirstChild;
[general_param]=load_var_general(childNode);

%% get param associated to a given embryo

params=[];
allListItems = docNode.getElementsByTagName('Embryo');
for k = 0:allListItems.getLength-1
    embIdx = k + 1;
    thisListItem = allListItems.item(k);
    childNode = thisListItem.getFirstChild;
    params{embIdx}=load_var_sub(childNode);
end

%% save the structure generated as template

[p,n]=fileparts(filexml);
mf=fullfile(p,[n '.mat']);
saveVarsMat_new.general_param = general_param;
saveVarsMat_new.params = params;
save(mf,'saveVarsMat_new');

clear all