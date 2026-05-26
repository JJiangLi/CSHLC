%-----------------------------------变量说明-------------------------------%
% I            原图像的double形式
% dW           方向模板宽度
% sW           平滑模板宽度
% N            平滑次数
% dGauMeanI    第N次方向高斯滤波后的结果  
% ctrEdgMap    经过膨胀处理后的边缘结果
%-------------------------------------------------------------------------%

function [dGauMeanI,ctrEdgMap] = edgeRegionSM(I,dW,sW,N)

    dirNum=8;                                      %模板方向数量
    [dKerList]=getDirDetKernel(dW,dirNum);         %返回方向模板列表
    [sKerList]=getSearchDir(sW,dirNum);            %返回平滑模板列表

    [row,col]=size(I);
%   NdirMap=zeros(row,col,N);    
    NedgMap=zeros(row,col,N);
    
    for m=1:N
        [dirMap]=detecDir(I,dKerList);             %返回每个像素点对应的最大响应的方向
        [I]= calDirMean(I,dirMap,sKerList,sW);     %返回方向高斯滤波的结果
        edgI=edge(I,'canny');                      %边缘提取，返回值0与1
%       NdirMap(:,:,m)=dirMap; 
        NedgMap(:,:,m)=edgI;                       %存放每次方向高斯平滑后提取的边缘
    end
    dGauMeanI=I;                                   %第N次方向高斯滤波后的结果
    
    %膨胀
    edgSelMap=ones(row,col);
    kernel=ones(3,3);
    for k=1:N-1
        edgIpd = imfilter(NedgMap(:,:,k),kernel,'replicate');  
        edgSelMap=edgSelMap&edgIpd;
    end

    ctrEdgMap=double(edgI&edgSelMap); 
    ctrEdgMap(ctrEdgMap~=0)=Inf;
    
%   [showEdgMap]= showConvKernel(dKerList,dW,dWNumber);  %3种形状
%   [showEdgMap]= showConvKernel(sKerList,sW,sWNumber);
    
end



%% 方向高斯滤波
function [DirGauMap]= calDirMean(I,dirMap,sKerList,sW)

    [row,col]=size(I);   
    [~,~,len]=size(sKerList);
    
    DirGauMap=zeros(row,col);
    ConvMap=zeros(row,col,len);
    gausFilter=fspecial('gaussian',[sW, sW],1.5);
    
    for n=1:len
        kernel=sKerList(:,:,n).*gausFilter;  
        kernel=kernel/sum(kernel(:));
        ConvMap(:,:,n)=imfilter(I,kernel, 'replicate'); 
    end  
    
    for i=1:row
        for j=1:col
            sel=dirMap(i,j);
            DirGauMap(i,j)=ConvMap(i,j,sel);
        end
    end

end


%% 方向检测
function [dirMap]=detecDir(I,decKerList)

    [row,col]=size(I);    
    [~,~,len]=size(decKerList); 
    ConvMap=zeros(row,col,len);
    
    for n=1:len
        kernel=decKerList(:,:,n);
        w=sum(abs(kernel(:)))/2;
        kernel=kernel/w;
        ConvMap(:,:,n)=imfilter(I,kernel, 'replicate'); 
    end        

    dirMap =zeros(row,col);   
    for x=1:row
        for y=1:col
            Dir=-1;
            Val=0;
            for k=1:len
                P=abs(ConvMap(x,y,k));
                if P>=Val
                   Val=P;
                   Dir=k;
                end
            end
            dirMap(x,y)=Dir;           
        end        
    end
    
end

 %% 方向均值点列表                          
function [sKerList]=getSearchDir(sW,dirNum)

    sR=(sW-1)/2;
    sKerList=zeros(sW,sW,dirNum); 
    edgMapTemp=zeros(sW*2+1,sW*2+1);
    edgMapTemp(:,sW+1)=2;  
    
    for k=0:dirNum
        theta=k*(180/dirNum); 
        %需要从大的块中间切小的模板
        edgMap=imrotate(edgMapTemp, theta-45,'nearest','crop');
        sKerList(:,:,k+1)=edgMap(sW+1-sR:sW+1+sR,sW+1-sR:sW+1+sR); 
    end
    
end

%% 方向检测模板
function [decKerList]=getDirDetKernel(dW,dirNum)

    dR=(dW-1)/2;                           %半径
    decKerList=zeros(dW,dW,dirNum*3);      %方向模板列表
    
    edgMapTemp=zeros(dW*2+1,dW*2+1);       %从大的中剪小的，防止形变
    %设置第一个方向模板
    edgMapTemp(1:dW+1,1:dW+1)=1; 
    edgMapTemp(dW+1:end,dW+1:end)=-1;
    edgMapTemp(dW+1,dW+1)=0;
    
    for k=0:dirNum-1
        theta=k*(180/dirNum);               %旋转角度
        %需要从大的块中间切小的模板
        edgMap=imrotate(edgMapTemp, theta,'nearest','crop');
        decKerList(:,:,k+1)=edgMap(dW+1-dR:dW+1+dR,dW+1-dR:dW+1+dR); 
    end    
    
end


%% 显示模板的内容
function [showEdgMap]= showConvKernel(kernelList,W,dirNumber)
    
    numW=dirNumber;                        %总模板数量
    colN=8;                                %列数量
    rowN=ceil(numW/colN);                  %行数量
    
    addW=2;
    maxW=W+addW;
    showEdgMap=zeros(maxW*rowN,maxW*colN); %显示图

    EdgTempList=zeros(maxW,maxW,numW);
    
    pw=(maxW-W)/2;
    for k=1:dirNumber
        edgCut=kernelList(:,:,k);
        EdgTempList(:,:,k)=padarray(edgCut,[pw,pw],200);
    end
    
    cnt=0;
    for i=0:rowN-1
        for j=0:colN-1
            cnt=cnt+1;
            showEdgMap(i*maxW+1:i*maxW+maxW,j*maxW+1:j*maxW+maxW)=(EdgTempList(:,:,cnt)+1)*100;
        end
    end
    
    edgMap=kron(showEdgMap,ones(5,5));       
    figure();imshow(uint8(edgMap));

end