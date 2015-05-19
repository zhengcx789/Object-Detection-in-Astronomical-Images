function outLabel=deblend_Multi(imlocal,imlocal_BW,DEBLEND_MINCONT,DEBLEND_NTHRESH)
% This function is for deblending the touching or overlapping objects by
% using the multi-threshold method,and using the watershed method to assign 
% the outlying pixels which have a flux lower than the separation thresholds 
% to their own proper objects
%
% Input: imlocal = the local smoothed image (It is a local area including 
%                  a object which needs to be checked and applied the deblending processing)
%        imlocal_BW = the local binary image(local masks) correspoinding to the same region of "imlocal",
%                     the value 1 in it means that this pixel belongs to an object) 
%        DEBLEND_MINCONT = the number of deblending sub-thresholds
%        DEBLEND_NTHRESH = the minimum contrast parameter for deblending
% Output:
%          outLabel = the local label matrix of deblended objects,different values 
%                     large than 0 in this matrix represent different objects                    

objects=imlocal.*imlocal_BW;
maxarea=[];
[rIndex,cIndex]=find(imlocal_BW==1);
thresh0=inf;
for i1=1:size(rIndex,1)
    if  objects(rIndex(i1,1),cIndex(i1,1))<thresh0
        thresh0=objects(rIndex(i1,1),cIndex(i1,1));
    end
end
peak=max(objects(:));
p=1;

%%%%%%%%%%%%%%%% Build tree(from bottom to top)%%%%%%%%%%%%%%%%
for k=1:DEBLEND_NTHRESH
    multithresh(p)=thresh0*((peak/thresh0)^(k/DEBLEND_NTHRESH)); % The levels are spaced exponentially.
    branch=(objects>=multithresh(p));
    tree{p}=branch;   % The larger p demonstrates the higher branch of the tree
    p=p+1;
end

%%%%%%%%%%%%%%%% Find delend branch(from top to bottom) %%%%%%%%%%%%%%%%
deblendtag=zeros(size(tree));
for j=size(tree,2):-1:1
    [branchlabel,branchnum]=bwlabel(tree{j},8);
    if branchnum>=2
        for i=1:branchnum
            smallmask=(branchlabel==i);
            smallobj=smallmask.*objects;
            intens(i)=sum(smallobj(:));
        end
        intensratio=intens./sum(objects(:));
        inds=find(intensratio>=DEBLEND_MINCONT);
        if size(inds,2)>=2
            deblendtag(j)=1;
            deblendlabel{j}=branchlabel; 
        end
    end
end

%%%%%%%%%%%%%%%% Find the deblended objects %%%%%%%%%%%%%%%%
La=find(deblendtag==1);
if isempty(La)
    outLabel=imlocal_BW;
elseif size(La,2)==1
    maxarea=deblendlabel{La(1)};  
elseif size(La,2)>1
    Laend=size(La,2);
    localmax{Laend}=deblendlabel{La(end)}; 
    for m=Laend-1:-1:1
        localmax{m}=deblendlabel{La(m)};
        inter=(localmax{m}>0).*(localmax{m+1}>0); 
        interind=find(inter==1);
        if ~isempty(interind) 
            removeL=unique(localmax{m}(interind)); 
            [localmaxM1,Lnum1]=bwlabel(localmax{m+1},8); 
            if Lnum1>size(removeL,1) 
                for rl=1:size(removeL,1)
                    localmax{m}(localmax{m}==removeL(rl))=0;
                end
                localmax{m}=localmax{m}+localmax{m+1};
            end
            
        else
            localmax{m}=localmax{m}+localmax{m+1}; 
            
        end
    end
    maxarea=localmax{m};
end

if ~isempty(maxarea)
    [ma_L,ma_num]=bwlabel(maxarea,8);
    for ii=1:ma_num
        smallmask=(ma_L==ii);
        smallobj=smallmask.*objects;
        intens(ii)=sum(smallobj(:));
    end
    intensratio=intens./sum(objects(:));
    ma_inds=find(intensratio<DEBLEND_MINCONT);
    if ~isempty(ma_inds)
        for remove_i=1:size(ma_inds,2)
            ma_L(ma_L==ma_inds(remove_i))=0; 
        end
    end
    maxlocal=(ma_L>0);  % the debelended objects
     
    %%%%%%%%%%%%%%  Assign the outlying pixels which have a flux lower than
    %%%%%%%%%%%%%%  the separation thresholds to their own proper objects
    [la,numm]=bwlabel(maxlocal,8);
    growmaxlocal=imdilate(la,strel('disk',2));
    outLabel=double(growmaxlocal).*double(imlocal_BW);
    if sum(objects(:))>2000
        D1=bwdist(maxlocal);
        D=imimposemin(D1,maxlocal);
        L1=watershed(D);
        outLabel=double(L1).*double(imlocal_BW);
    end
end


