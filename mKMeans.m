function [resMap] = mKMeans(I,C)
    
    [Idx,cent]=kmeans(single(I(:)),C);
    CMap=reshape(Idx,size(I));
    [~,ind]=sort(cent);
    [~,ind3]=sort(ind);
    [row,col]=size(I);
    
    resMap=zeros(row,col);
    for i=1:row
        for j=1:col
            nexLabel=ind3(CMap(i,j));
            resMap(i,j)=nexLabel;
        end
    end

end