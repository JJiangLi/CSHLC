%--------------------------------变量说明----------------------------------%
% imgName  待分割图像
% CT       待分割图像的地面真值
% C        待分割图像的类数
%-------------------------------------------------------------------------%

clear; close all; clc;

imageName={};
imageGT={};
C=[];

for k=1:length(imageName)
    datestr(now)                              %以字符串的形式返回当前时间
    [Acc]=autoMain(imageName{k},imageGT{k},C(k));
end


function [Acc]=autoMain(imgName,GT,C)

    disp(['正在跑数据：',num2str(imgName)]);
    I=imread(['img/',imgName,'.bmp']);
    GT=imread(['img/',GT,'.bmp']);
    saveName=imgName;
    saveExd='.png';
    savePath='res/';
    
    tic
    [finalI,saveParam]=CSHLC(I,C);
    tim=toc;
    
    [grayMap,rgbMap,Acc,Kappa] = calAccKappa(finalI,GT);
    %保存图像和准确率，时间
    saveAcc=num2str(round(Acc*1000));
    saveTim=num2str(round(tim*1000));
    saveKappa=num2str(round(Kappa*1000));
    saveGray=[savePath,'gray/',saveName,'_',saveAcc,'A_',saveKappa,'K_',saveTim,'ms',saveExd];
    saveRGB=[savePath,'rgb/',saveName,'_',saveAcc,'A_',saveKappa,'K_',saveTim,'ms',saveExd];
    imwrite(grayMap,saveGray);
    imwrite(rgbMap,saveRGB);
    
    file=fopen([savePath,'parameter.txt'],'a+');                                                      
    fprintf(file,' %s',saveRGB);
    fprintf(file,' %s \n ',saveParam);
    fclose(file);
    
    disp(['保存图片和参数成功,准确率：',num2str(Acc),' Kappa: ',num2str(Kappa) ,' 时间：',num2str(tim)]);

end

function [grayMap,rgbMap,Acc,Kappa] = calAccKappa(resMap,groundT)

    % 颜色已配对的图计算准确率,保存图片
    [row,col]=size(resMap);
    [valList,~]=sort(unique(groundT(:)));
 
    grayMap=zeros(row,col,'uint8');               %生成灰度图  
    for i=1:row
        for j=1:col
            grayMap(i,j)=valList(resMap(i,j));
        end
    end
    
    color=[ 61, 38,168; %深蓝
            39,150,235; %浅蓝
           128,203, 88; %绿            
           249,250, 20; %黄
           208,191, 39; %棕
           205,  0,  0; %红
        ];
    rgbMap=zeros(row,col,3,'uint8');
    for i=1:row
        for j=1:col
            cur=resMap(i,j);
            rgbMap(i,j,1)=color(cur,1);
            rgbMap(i,j,2)=color(cur,2);
            rgbMap(i,j,3)=color(cur,3);
        end
    end
    
    pe =0;
    for j =1:length(valList)
        tmp_1=(grayMap==valList(j));
        tmp_2=(groundT==valList(j));
        pe = pe + sum(tmp_1(:))*sum(tmp_2(:));
    end
    total=row*col;
    pe = pe/(total*total);
    
    accMap = (grayMap==groundT);
    acc = sum(accMap(:))/total;
    Kappa = (acc-pe)/(1-pe);
    
    Acc=100*acc; %计算准确率 
    
%     Acc=100*sum(sum(grayMap==groundT))/(row*col); %计算准确率 
 
    
end

%% 原算法核心    
function [finalI,saveParam]=CSHLC(I,C)

    I=double(I);
    
    % 边缘提取
    dW=7;       %方向模板宽度7
    sW=dW;      %平滑模板宽度7
    N=6;        %平滑次数6
    [dGauMeanI,ctrEdgMap] = edgeRegionSM(I,dW,sW,N);
    % dGauMeanI    第N次方向高斯滤波后的结果  
    % ctrEdgMap    经过膨胀处理后的边缘结果

    % 边缘约束平滑
    ctrW=9;        %滑动窗口大小
    [smI] = ctrSmooth(dGauMeanI,ctrEdgMap,C,ctrW);
    
    % CRF 
    CRF_iter=4;
    CLabelMap=mKMeans(smI,C);
    [CRFLabelMap] =getPixelCRF(I,CLabelMap,C,CRF_iter );
    % CRFLabelMap 返回的是使用MRF重新标记的结果

    W=13;
    minNum=1400;
    [corrMap]=selMaxCorr(I,CRFLabelMap,CLabelMap);
    [corrMap2]=regionGrowCorr(I,corrMap,ctrEdgMap,W,C,minNum);
 
    [corrMap2] =getPixelCRF(I,corrMap2,C,CRF_iter);
    [finalI]=selMaxCorr(I,corrMap2,CLabelMap);

    % 保存运行参数
    saveParam=[' dW=',num2str(dW),'  sW=',num2str(sW),'  N=',num2str(N)];
    saveParam=[saveParam,' ctrW=',num2str(ctrW),'  CRFiter=',num2str(CRF_iter)];
    saveParam=[saveParam,' C=',num2str(C),'  W=',num2str(W),'  minNum=',num2str(minNum)];
end


