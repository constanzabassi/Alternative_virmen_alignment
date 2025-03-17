function allcells_st = make_allcells (resp_tr)
allcells_st =[];
for d = 1:length(resp_tr)
    allcells = {};
    for c = 1:size(resp_tr{1,d},2)
        allcells(c).trials = squeeze(resp_tr{1,d}(:,c,:));
    end
    allcells_st= [allcells_st,allcells];
end