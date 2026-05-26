%--------------------------------变量说明----------------------------------%
% I           待分割图像
% corrMap     像素群计数对比的修正结果 
% selEdgMap   经过膨胀处理后的边缘结果
% W           生长的最大窗口宽度
%-------------------------------------------------------------------------%
function [corrMap2]=regionGrowCorr(I,corrMap,selEdgMap,W,C,minNum)

    edgI=edge(corrMap,'canny');       %对角点修复后的类标图进行边缘提取
    edgI2=edgI&selEdgMap;             %二次筛选
    
    [reLabelMap_CRF,labelCnt_CRF]= getReLabelMap(corrMap);
    [LabelNumMap_CRF,LabelMeanMap_CRF]= getNumMeanMap(I,reLabelMap_CRF,labelCnt_CRF);
    
    edgCLabelMap=corrMap;
    edgCLabelMap(edgI2~=0)=Inf;       %边界标记
    
    selMap=(LabelNumMap_CRF<minNum);  %统计像素群像素个数小于指定个数的像素  
    
    [row,col]=size(corrMap);
    fr=(W-1)/2;
    
    for i=1:row
        for j=1:col
            if edgCLabelMap(i,j)~=Inf && selMap(i,j)
                X0= max(1,i-fr);
                Xn= min(row,i+fr);
                Y0= max(1,j-fr);
                Yn= min(col,j+fr);
                srchCut=edgCLabelMap(X0:Xn,Y0:Yn);         %
                MeanCut=LabelMeanMap_CRF(X0:Xn,Y0:Yn);
                [mostLabel]=findMostSimilarLabel(srchCut,MeanCut,C,i-X0+1,j-Y0+1);
                edgCLabelMap(i,j)=mostLabel;
            end
        end
    end
    
    corrMap2=edgCLabelMap;
    edgI3=edge(corrMap2,'canny');
    
    corrMap2(edgI2~=0)=0;
    corrMap(edgI2==0)=0;

    corrMap2=corrMap2+corrMap;


end

function [mostLabel]=findMostSimilarLabel(cut,MeanCut,C,i,j)

    adj=[-1,0;0,-1;1,0;0,1];
    
    [row,col]=size(cut);
    labelVec=zeros(C,1);   
    sumVec=zeros(C,1);
    roadMap=zeros(row,col);
    
    Xarray=zeros(row*col,1); 
    Yarray=zeros(row*col,1);
    AftPoint=1;                  
    BefPoint=1;                  
    Xarray(AftPoint)=i;
    Yarray(AftPoint)=j;

    while BefPoint <= AftPoint
        x=Xarray(BefPoint);
        y=Yarray(BefPoint);
        BefPoint=BefPoint+1;
        roadMap(x,y)=1;
        labelVec(cut(x,y),1)=labelVec(cut(x,y),1)+1;  
        sumVec(cut(x,y),1)=sumVec(cut(x,y),1)+MeanCut(x,y);
        for k=1:length(adj)
           xx=x+adj(k,1);
           yy=y+adj(k,2);
           if xx>=1 && xx<=row && yy>=1 && yy <=col && cut(xx,yy)~=Inf
               if (roadMap(xx,yy)==0)  
                   AftPoint=AftPoint+1;
                   Xarray(AftPoint)=xx;
                   Yarray(AftPoint)=yy;
                   roadMap(xx,yy)=2;                                   
               end
           end
        end
    end
    
    [maxnum,mostLabel]=max(labelVec);
    
    MeanVec=sumVec./labelVec;
    MeanVec(cut(i,j),1)=Inf;
    
    disVec=abs(MeanVec-MeanCut(i,j));
    [minnum,L2]=min(disVec);
    
    if abs(mostLabel-cut(i,j))>1
        mostLabel=L2;
    end

end

%% 统计每个像素对应的类中像素个数
function [LabelNumMap,LabelMeanMap]=getNumMeanMap(I,spLabelMap,spLabelNum)

    [row,col]=size(spLabelMap);
    LabelNumVec=zeros(spLabelNum,1);
    LabelSumVec=zeros(spLabelNum,1);
    for i=1:row
        for j=1:col
            label=spLabelMap(i,j);
            if label~=0 && label~=Inf
                LabelNumVec(label)=LabelNumVec(label)+1;
                LabelSumVec(label)=LabelNumVec(label)+1;
            end
        end
    end
    LabelMeanVec=LabelSumVec./LabelNumVec;
    
    LabelNumMap=zeros(row,col);
    LabelMeanMap=zeros(row,col);
    for i=1:row 
        for j=1:col
            tempLabel=spLabelMap(i,j);
            if tempLabel~=0 && tempLabel~=Inf
                LabelNumMap(i,j)=LabelNumVec(tempLabel);
                LabelMeanMap(i,j)=LabelMeanVec(tempLabel);
            end
        end
    end
   
end


%% 像素群类标统计
function [reEdgLabelMap,labelCnt]= getReLabelMap(edgCLabelMap)

    adj=[-1,0;0,-1;1,0;0,1];    
    
    [row,col]=size(edgCLabelMap);
    globalRoadMap=zeros(row,col);
    reEdgLabelMap=zeros(row,col)+Inf;
    labelCnt=0;    
    Xarray=zeros(row*col,1); 
    Yarray=zeros(row*col,1);

    for i=1:row 
        for j=1:col
            if  edgCLabelMap(i,j)~=Inf && globalRoadMap(i,j)==0 
                globalRoadMap(i,j)=1;
                labelCnt=labelCnt+1;                
                AftPoint=1;                   
                BefPoint=1;                   
                Xarray(AftPoint)=i;
                Yarray(AftPoint)=j;                
                while BefPoint <= AftPoint
                    x=Xarray(BefPoint);
                    y=Yarray(BefPoint);
                    BefPoint=BefPoint+1;
                    
                    currentLabel=edgCLabelMap(x,y);
                    globalRoadMap(x,y)=1;
                    reEdgLabelMap(x,y)=labelCnt;                    
                    for k=1:length(adj)
                       xx=x+adj(k,1);
                       yy=y+adj(k,2);
                       if xx>=1 && xx<=row && yy>=1 &&yy <=col && edgCLabelMap(xx,yy)~=Inf
                           if globalRoadMap(xx,yy)==0 && edgCLabelMap(xx,yy)==currentLabel 
                               AftPoint=AftPoint+1;
                               Xarray(AftPoint)=xx;
                               Yarray(AftPoint)=yy;
                               globalRoadMap(xx,yy)=2;                                   
                           end
                       end
                    end
                end
            end
        end
    end
    
end