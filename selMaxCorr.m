%--------------------------------变量说明----------------------------------%
% I              待分割图像
% CRFLabelMap    CRF修复后的类标图
% CLabelMap      K_Means聚类之后的类标图
%-------------------------------------------------------------------------%
function [corrMap]=selMaxCorr(I,CRFLabelMap,CLabelMap)

    [reLabelMap_CRF,labelCnt_CRF]= getReLabelMap(CRFLabelMap);
    [LabelNumMap_CRF,~]= getNumMeanMap(I,reLabelMap_CRF,labelCnt_CRF);
   
    [reLabelMap_C,labelCnt_C]= getReLabelMap(CLabelMap);
    [LabelNumMap_C,~]= getNumMeanMap(I,reLabelMap_C,labelCnt_C);
    
    corrMap=CRFLabelMap;
    minNum=10;
    [row,col]=size(CRFLabelMap);
    for i=1:row
        for j=1:col
           if LabelNumMap_CRF(i,j)<minNum
              if  LabelNumMap_C(i,j)>LabelNumMap_CRF(i,j)
                  corrMap(i,j)=CLabelMap(i,j);
              end
           end
        end
    end

end


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
