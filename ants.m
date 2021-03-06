function [path,min_dist] = ants(p)


m=20;                                            %蚂蚁数量
Alpha=1;                                         %信息素重要程度参数
Beta=5;%5;                                          %启发式重要程度参数
Rho=0.3;                                         %信息素蒸发参数
G=400;                                           %最大迭代次数
Q=100;%100;                                           %信息素增加强度系数
C=p;                                             %需要规划的点的坐标

%第一步：变量初始化
n=size(C,1);                                     %问题的规模，即城市数目
D=zeros(n,n);                                    %任何两个城市距离间隔矩阵
for i=1:n
    for j=1:n
        if i~=j
            D(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
        else
            D(i,j)=eps;
        end
        D(j,i)=D(i,j);
    end
end
Eta=1./D;                                     %Eta为启发因子，设为距离的倒数
Tau=ones(n,n);                                %信息素矩阵
Tabu=zeros(m,n);                              %信息素矩阵
NC=1;                                         %迭代计数器
R_Best=zeros(G,n);                            %各代最佳路线
L_Best=inf.*ones(G,1);                        %各代最佳路线长度,inf正无穷
figure(1);                                    %优化解
while NC<=G
    %第二步：将m只蚂蚁放到n个城市上
    Randpos=[];
    for i=1:(ceil(m/n))
        Randpos=[Randpos,randperm(n)];
    end
    Tabu(:,1)=(Randpos(1,1:m))';
    %第三步：m只蚂蚁按照概率函数选择下一座城市，完成各自周游
    for j=2:n
        for i=1:m
            visited=Tabu(i,1:(j-1));                      %已访问的城市
            J=zeros(1,(n-j+1));                           %待访问的城市
            P=J;                                          %待访问城市的选择概率分布
            Jc=1;
            for k=1:n
                if isempty(find(visited==k,1))
                    J(Jc)=k;
                    Jc=Jc+1;
                end
            end
            %计算待选城市概率分布
            for k=1:length(J)
                P(k)=(Tau(visited(end),J(k))^Alpha*(Eta(visited(end),J(k))^Beta));
            end
            P=P/(sum(P));
            %按概率选取下一个城市
            Pcum=cumsum(P);
            Select=find(Pcum>=rand);
            to_visit=J(Select(1));
            Tabu(i,j)=to_visit;
        end
    end
    
    if NC>2
        Tabu(1,:)=R_Best(NC-1,:);%这儿看似是一只蚂蚁每一次跑的结果都被舍弃了，但是实际上是，之前的最优结果能一直对后续蚂蚁产生信息素增强效果
    end
    L=zeros(m,1);
    for i=1:m
        R=Tabu(i,:);
        for j=1:(n-1)
            L(i)=L(i)+D(R(j),R(j+1));
        end
        L(i)=L(i)+D(R(1),R(n));
    end
    L_Best(NC)=min(L);
    pos=find(L==L_Best(NC));
    R_Best(NC,:)=Tabu(pos(1),:);
    %第五步：更新信息素
    Delta_Tau=zeros(n,n);
    for i=1:m
        for j=1:(n-1)
            Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
    end
    Tau=(1-Rho).*Tau+Delta_Tau;
    %第六步：禁忌表清零
    Tabu=zeros(m,n);
    %历代最优路线
    str = ['第' num2str(NC) ' 次迭代最优距离: ' num2str(L_Best(NC))];
    disp(str);
    for i=1:(n-1)
        plot([C(R_Best(NC,i),1),C(R_Best(NC,i+1),1)],[C(R_Best(NC,i),2),C(R_Best(NC,i+1),2)],'bo-');
        hold on
    end
    plot([C(R_Best(NC,n),1),C(R_Best(NC,1),1)],[C(R_Best(NC,n),2),C(R_Best(NC,1),2)],'ro-');
    %title(['优化最优距离：',num2str(L_Best(NC))]);
    title(str);
    hold off;
    pause(0.005);
    NC=NC+1;
end
%第七步：输出结果
pos=find(L_Best==min(L_Best));
Shortest_Route=R_Best(pos(1),:);                             %最佳路线
Shortest_Length=L_Best(pos(1));                              %最佳路线长度
figure(2),
plot(L_Best)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')