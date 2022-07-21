function I=ea_gen_pseudoimprovements(list,sweetspots,weights)

if ~exist('weights','var')
    weights=ones(length(sweetspots),1);
end


% gather sweetspot coords:
val=zeros(length(list),length(sweetspots));
for sws=1:length(sweetspots)
    [pth,fn,ext]=fileparts(sweetspots{sws});
    switch ext
        case {'.nii','.gz'} % roi / nucleus
            roi=ea_load_nii(sweetspots{sws});
            [xx,yy,zz]=ind2sub(size(roi.img),find(roi.img(:)));
            XYZ=[xx,yy,zz,ones(size(xx,1),1)]';
            roiMm=roi.mat*XYZ; % port to mm (world) space
            roiMm=roiMm(1:3,:)';
            roiVals=roi.img(find(roi.img(:)));
        case {'.mat'} % fibertract
            roi=load(sweetspots{sws});
            roiMm=double(roi.fibers(:,1:3));
            roiVals=ones(size(roiMm,1),1);
            roi.voxsize=mean(abs(diff(roiMm)),1); % fake a rough spatial resolution
    end

    KD=KDTreeSearcher(roiMm);

    ea_dispercent(0,['Calculating scores for sweetspot #',num2str(sws)])
    for vta=1:length(list)
        stim=ea_load_nii(list{vta});
        [xx,yy,zz]=ind2sub(size(stim.img),find(stim.img(:)>0));
        XYZ=[xx,yy,zz,ones(size(xx,1),1)]';
        stimMm=stim.mat*XYZ; % port to mm (world) space
        stimMm=stimMm(1:3,:)';
        stimVals=stim.img(find(stim.img(:)>0));


        [IDX,D]=knnsearch(KD,stimMm);
        connidx=D<(mean(roi.voxsize)*1.5);
        if ~any(connidx)
            val(vta,sws)=nan;
        else
            try
            val(vta,sws)=ea_nansum(roiVals(IDX(connidx)).*stimVals(connidx))./...
                ea_nansum(stimVals(~connidx));
            catch
                            val(vta,sws)=nan;
            end
        end
        
        ea_dispercent(vta/length(list));
    end
    ea_dispercent(1,'end');
end

I=val.*repmat(weights,size(val,1),1);
I=ea_nanmean(I,2);






