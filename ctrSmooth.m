%--------------------------------变量说明----------------------------------%
% dirGauMeanI  第N次方向高斯滤波后的结果
% edgI         经过膨胀处理后的边缘结果
% C            待分割图像的类数
% W            区域生长窗口大小
%-------------------------------------------------------------------------%
function [finalI] = ctrSmooth(dirGauMeanI,edgI,C,W)
        
    fr=(W-1)/2;                                %半径
    [row,col]=size(dirGauMeanI);
    
    sigma=2;
    GauTemp=fspecial('gaussian',[W, W],sigma); %获取高斯滤波核 
    
    smoothMap=zeros(row,col);
    for i=1:row
        for j=1:col
            if edgI(i,j)==0                    %如果该点不是边缘
                X0= max(1,i-fr);
                Xn= min(row,i+fr);
                Y0= max(1,j-fr);
                Yn= min(col,j+fr);
                edgCut=edgI(X0:Xn,Y0:Yn);      %“截取”边缘提取图像
                ICut=dirGauMeanI(X0:Xn,Y0:Yn); %“截取”方向高斯滤波的图像
                %截取高斯滤波算子
                tempCut=GauTemp(X0-i+1+fr:Xn-i+1+fr,Y0-j+1+fr:Yn-j+1+fr);      
                %实质是匀质区域滤波
                [smoothVal]=smoothUsingEdg(ICut,edgCut,tempCut, C,i-X0+1,j-Y0+1); 
                smoothMap(i,j)=smoothVal;
            end
        end
    end
    
    smoothMap(smoothMap==0)=Inf;               %恢复边缘点
    [finalI]=transEdgToLabel(smoothMap,dirGauMeanI);

end

function [resMap]=transEdgToLabel(edgLabelMap,Img)

    adj=[-1,0;0,-1;1,0;0,1];
    [row,col]=size(edgLabelMap);
    resMap=edgLabelMap;

    [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap);
    while(hasEdg)                             %当存在边缘时
        for i=1:row
            for j=1:col
                if outerEdgMap(i,j)==1        %当该点为边缘时
                   mindis=Inf;
                   mink=1;
                   for k=1:length(adj)
                        xx=i+adj(k,1);
                        yy=j+adj(k,2);
                        if xx>=1 && xx<=row && yy>=1 && yy <=col && edgLabelMap(xx,yy)~=Inf
                            d=abs(Img(xx,yy)-Img(i,j));
                            if d<mindis
                               mindis=d;
                               mink=k;
                            end
                        end
                   end
                   resMap(i,j)=edgLabelMap(i+adj(mink,1),j+adj(mink,2));
                end
            end
        end  
        edgLabelMap=resMap;
        [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap);  
    end 
end

function [outerEdgMap,hasEdg]=findOuterEdg(edgLabelMap)

    adj=[-1,0;0,-1;1,0;0,1;1,1;1,-1;-1,1;-1,-1];
    [row,col]=size(edgLabelMap);
    outerEdgMap=zeros(row,col);
    
    hasEdg=0;                                %标志位
    for i=1:row
        for j=1:col
            if edgLabelMap(i,j)==Inf         %如果该点是边缘，则对其8邻域进行搜索
                for k=1:length(adj)
                   xx=i+adj(k,1);
                   yy=j+adj(k,2);
                   if xx>=1 && xx<=row && yy>=1 && yy <=col && edgLabelMap(xx,yy)~=Inf
                      outerEdgMap(i,j)=1;    %如果该点邻域内没有其他边缘，则该点值为1
                      if hasEdg==0
                         hasEdg=1; 
                         break;
                      end
                   end
                end
            end
        end
    end

end
    
function [smoothVal]=smoothUsingEdg(ICut,edgCut,GauTemp,C,i,j)

    adj=[-1,0;0,-1;1,0;0,1];  %4邻域
    [row,col]=size(edgCut);   
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
        for k=1:length(adj)
           xx=x+adj(k,1);
           yy=y+adj(k,2);
           if xx>=1 && xx<=row && yy>=1 && yy <=col && edgCut(xx,yy)~=Inf
               if (roadMap(xx,yy)==0) 
                   AftPoint=AftPoint+1;
                   Xarray(AftPoint)=xx;
                   Yarray(AftPoint)=yy;
                   roadMap(xx,yy)=2;                                   
               end
           end
        end
    end
    GauTemp(roadMap==0)=0;    %像素点是边缘则对应的高斯滤波位置值为0，即不滤波
    res=(ICut.*GauTemp);
    smoothVal = sum(res(:))/sum(GauTemp(:));

end
