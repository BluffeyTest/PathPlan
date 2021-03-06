function [path,min_dist] = findBestPath(p)
% [path,minDist] = findBestPath(Points)
% path 路径点，按顺序排列
% min_dist 最短回环路径距离之和
% points 路径点的集合，按纵列排列，第一列是x，第二列是y

[m,n]=size(p);
if n ~= 2
    p = p';
end

%距离矩阵
distance = getDistanceMatrix(p);
OneEdgeTopBoder = max(max(distance));    % 最小单边距离
PathSum = OneEdgeTopBoder*m;            %记录的当前的路径总长度

% 迭代500次，每次从固定点出发,暂且就是第一个点吧
figure(1);
[path,min_dist,curOneEdgeTopBoder]=getPath(distance,OneEdgeTopBoder)
OneEdgeTopBoder = curOneEdgeTopBoder;
for i=2:50000
    
    [curpath,dist,curOneEdgeTopBoder]=getPath(distance,OneEdgeTopBoder);
    if dist<min_dist && dist~=-1
        min_dist = dist
        path = curpath;
        OneEdgeTopBoder = curOneEdgeTopBoder;
    end

    %在这儿也绘制一把
    path_points = zeros(m,n);
    for i=1:m
        path_points(i,1) = p(path(1,i),1);
        path_points(i,2) = p(path(1,i),2);
    end
    plot(path_points(:,1),path_points(:,2),'o-');
    hold on;
    plot([path_points(1,1),path_points(m,1)],[path_points(1,2),path_points(m,2)],'ro-');
    title(['优化最优距离：',num2str(min_dist)]);
    hold off;
    pause(0.005);
end
    


end



%% 单次迭代
function [path,dist,curOneEdgeTopBoder]=getPath(distance,OneEdgeTopBoder)
    b = 0;%是否没法完成计算需要退出
    [m,~] = size(distance);
    n = 54;
    curPath =  zeros(1,m);%当前顺序路径点在路径矩阵中的顺序
    curPath(1,1) = n;
    curUse = zeros(1,m);    %当前路径点当中哪个位置的点已经被使用了
    curUse(1,n) = 1;
    
    for j=2:m
        %寻找下一个点的基本思路是，当下一个点到当前点的距离不大于单边上界时，按1-单边距离/上界距离作为基础概率选择基础概率最大的点
        %考虑局部随机贪心，比如从当前最近的5个点当中随机选择一个
        prob=zeros(1,m);
        for i = 1:m
            if (curUse(1,i)==0)  && (i~=n) %&& (distance(n,i)<=OneEdgeTopBoder)
                prob(1,i) = 1-distance(n,i)/OneEdgeTopBoder;
                if prob(1,i) == 0
                    prob(1,i) = rand();%当为0的时候随机取一个概率
                end
                %prob(1,i) = rand() * prob(1,i);
            end
            
        end
        %贪心
        %n = find(prob == max(max(prob)));;

        %局部随机贪心
        prob_sort = sort(prob,2,'descend');
        num = m-sum(sum(curUse));
        if num>1
            num = 1;
        end
        r = rand(1,num);
        n1 = find(r==max(max(r)));
        n = find(prob == prob_sort(1,n1));


        [nm,nn] = size(n);
        if nm~=1 || nn~=1
            b = 1;
            dist=-1;
            break;
        end
        curUse(1,n) = 1;
        curPath(1,j) = n;
    end
    
    path = curPath;
    curOneEdgeTopBoder = OneEdgeTopBoder;
    if b==0
        curOneEdgeTopBoder = 0;
        dist = distance(curPath(1,1),curPath(1,m));
        for i=2:m
            d = distance(curPath(1,i-1),curPath(1,i));
            dist = dist+d;
            if d>curOneEdgeTopBoder
                curOneEdgeTopBoder = d;
            end
        end
        curOneEdgeTopBoder = curOneEdgeTopBoder*2;
        if curOneEdgeTopBoder>OneEdgeTopBoder
            curOneEdgeTopBoder = OneEdgeTopBoder;
        end
        
    end
    
end




%%
function [distance]=getDistanceMatrix(p)
%   [distance]=getDistanceMatrix(p)
%   获得距离矩阵
%   distance    距离矩阵
%   p   点集合，按纵列排列
[m,~] = size(p);
distance = zeros(m,m);
for i=1:m
    distance(i,i) = -1;
    for j=i+1:m
        distance(i,j) = sqrt((p(i,1)-p(j,1))^2 + (p(i,2)-p(j,2))^2);
        distance(j,i) = distance(i,j);
    end
end
end
