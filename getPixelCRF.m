%--------------------------------变量说明----------------------------------%
% labelMap  K_Means结果
% img       待分割图像
% maxiter   最大迭代次数
%-------------------------------------------------------------------------%
function [labelMap] =getPixelCRF(img,labelMap,C,maxiter )

    [row,col] =size(img);

    for iter=0:maxiter-1
        
        %这里采用的是像素点和3*3领域的标签相同与否来作为计算概率
        %收集上下左右斜等八个方向的标签
        label_u = imfilter(labelMap,[0,1,0; 0,0,0; 0,0,0],'replicate');
        label_d = imfilter(labelMap,[0,0,0; 0,0,0; 0,1,0],'replicate');
        label_l = imfilter(labelMap,[0,0,0;1,0,0;0,0,0],'replicate');
        label_r = imfilter(labelMap,[0,0,0;0,0,1;0,0,0],'replicate');
        label_ul = imfilter(labelMap,[1,0,0;0,0,0;0,0,0],'replicate');
        label_ur = imfilter(labelMap,[0,0,1;0,0,0;0,0,0],'replicate');
        label_dl = imfilter(labelMap,[0,0,0;0,0,0;1,0,0],'replicate');
        label_dr = imfilter(labelMap,[0,0,0;0,0,0;0,0,1],'replicate');
        p_c = zeros(C,row*col);
        
        %计算像素点8领域标签相对于每一类的相同个数
        for i = 1:C
            label_i = i * ones(row,col);
            temp = ~(label_i - label_u) + ~(label_i - label_d) + ...
                ~(label_i - label_l) + ~(label_i - label_r) + ...
                ~(label_i - label_ul) + ~(label_i - label_ur) + ...
                ~(label_i - label_dl) +~(label_i - label_dr);
            p_c(i,:) = temp(:)/8;%计算概率
        end
        p_c((p_c == 0)) = 0.001; %防止出现0
        
        %计算似然函数
        mu = zeros(1,C);         %均值
        sigma = zeros(1,C);      %方差
        %求出每一类的的高斯参数：均值方差
        for i = 1:C
            data_c = img((labelMap == i));
            mu(i) = mean(data_c);  %均值
            sigma(i) = var(data_c);%方差
        end
        
        p_sc = zeros(C,row*col);
        %计算每个像素点属于每一类的似然概率
        %为了加速运算，将循环改为矩阵一起操作
        for j = 1:C
            MU = repmat(mu(j),row*col,1);
            p_sc(j,:) = 1/sqrt(2*pi*sigma(j))*exp(-(img(:)-MU).^2/(2*sigma(j)));
        end
        %找到联合一起的最大概率最为标签，取对数防止值太小
        [~,labelMap] = max(log(p_c) + log(p_sc));
        labelMap = reshape(labelMap,[row,col]);
        
    end
    
    %[labelMap] = relabelByGrayASC(labelMap, img);
    [labelMap] = relabelByGrayASC(labelMap, img);
    

end


function [newLabelMap] = relabelByGrayASC(labelMap, I)
    
    [row,col]=size(I);
    
    numVec =unique(labelMap(:));
    len = length(numVec);
    sumVec =zeros(len,1);
    cntVec =zeros(len,1);
    for i=1:row
        for j=1:col
            label = labelMap(i,j);
            sumVec(label) = sumVec(label) + I(i,j);      %计算同一类像素灰度值的和
            cntVec(label) = cntVec(label)+1;             %计算同一类像素的个数
        end
    end
    sumVec = sumVec./cntVec;                             %计算同一类像素的平均灰度值

    [sortVec, indVec]=sort(sumVec);
    newLabelMap=zeros(row,col);
    for i=1:row
        for j=1:col
            label = labelMap(i,j);
            for k=1:len
                if label ==indVec(k)
                    newLabelMap(i,j)=k;
                    break;
                end
            end
        end
    end
       
%     figure();imagesc(newLabelMap)
%     figure();imagesc(labelMap)

end
