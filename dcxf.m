%约束矩阵:A
%目标函数系数向量：c
%约束右端向量：b
%目标函数取最小值时的自变量值：X
%目标函数的最小值：f

A=[1,-2,1,0;1,2,0,1];
b=[4 8];
c=[-1 -1 0 0];

format  rat     %可以让结果用分数输出
[m,n]=size(A);
E=1:m;E=E';              
F=n-m+1:n;F=F';
D=[E,F];             %创建一个一一映射，为了结果能够标准输出
X=zeros(1,n);        %初始化X
if n<m              %判断是否为标准型
    disp('不符合要求需引入松弛变量')
    flag=0;
else
    flag=1;
    B=A(:,n-m+1:n);   %找基矩阵
    cB=c(n-m+1:n);    %基矩阵对应目标值的c
   
   while flag 
     w=cB/B;         %计算单纯形乘子，cB/B=cB*inv(B),用cB/B的目的是，为了提高运行速度
     pbs=w*A-c;
     [z,k]=max(pbs);  % k作为进基变量下标 
     if(z<0.000000001)     
         flag=0;             %所有判别数都小于0时达到最优解
         disp('    已找到最优解');        
         xB=(B\b')';
         f=cB*xB';   
         for i=1:n
           mark=0;
             for j=1:m
                if (D(j,2)==i)
                 mark=1;
                 X(i)=xB(D(j,1));     %利用D找出xB与X之间的关系
                end
             end
             if mark==0
               X(i)=0;         %如果D中没有X(i),则X(i)为非基变量，所以X(i)＝0
             end
         end
          disp('基向量为:'); 
          disp(X);
          disp('目标函数值为:') ; 
          disp(f);
     else
         if(B\A(:,k)<=0)          % 如果B\A(;,k)中的每一个分量都小于零
               flag=0;
               disp(' 此问题不存在最优解');  %若B\A(:,k)的第k列均不大于0，则该问题不存在最优解
         else
               b1=B\b';
               temp=inf;
               for i=1:m
                   if ((A(i,k)>0) && (b1(i)/(A(i,k)+eps))<temp )
                       temp=b1(i)/A(i,k);   %找退基变量
                       r=i;
                   end                
               end
               B(:,r)=A(:,k);
               cB(r)=c(k); %确定换入换出变量后，相应的基矩阵及新基对应的目标值的c也相应改变
               D(r,2)=k;   %改变D中的映射关系
              
         end
     end
   end
end
