function rechuandaopde 
  %以下所用数据，除了t的范围我根据题目要求取到了20000，其余均从pdf中得来 
  a=0.00045;%a的取值 
  xspan=[0 10];%x的取值范围 
  tspan=[0 5400];%t的取值范围 
  ngrid=[1000 100];%分割的份数，前面的是t轴的，后面的是x轴的 
  f=@(x)37;%初值 
  g1=@(t)75;%边界条件一 
  g2=@(t)75;%边界条件二 
  [T,x,t]=pdesolution(a,f,g1,g2,xspan,tspan,ngrid);%计算所调用的函数 
  [x,t]=meshgrid(x,t); 
  mesh(x,t,T);%画图，并且把坐标轴名称改为x，t，T 
  xlabel('x') 
  ylabel('t') 
  zlabel('T') 
  T%输出温度矩阵 
  dt=tspan(2)/ngrid(1);%t步长 
  h3000=3000/dt; 

h9000=9000/dt; 
h15000=15000/dt;%3000,9000,15000下，温度分别在T矩阵的哪些行 
T3000=T(h3000,:) 
T9000=T(h9000,:) 
T15000=T(h15000,:)%输出三个时间下的温度分布 
%不再对三个时间下的温度-长度曲线画图，其图像就是三维图的截面 
 
%稳定性讨论,傅里叶级数法 
dx=xspan(2)/ngrid(2);%x步长 
sta=4*a*dt/(dx^2)*(sin(pi/2))^2; 
if sta>0,sta<2 
    fprintf('\n%s\n','有稳定性') 
else  
    fprintf('\n%s\n','没有稳定性') 
    error 
end 
 
%真实值计算 
[xe,te,Te]=truesolution(a,f,g1,g2,xspan,tspan,ngrid); 
[xe,te]=meshgrid(xe,te); 
mesh(xe,te,Te);%画图，并且把坐标轴名称改为xe，te，Te 
xlabel('xe') 
ylabel('te') 
zlabel('Te') 
Te%输出温度矩阵 
 
%误差计算 
jmax=1/dx+1;%网格点数 
[rms]=wuchajisuan(T,Te,jmax) 
rms%输出误差 

 
function [rms]=wuchajisuan(T,Te,jmax) 
for j=1:jmax 
    rms=((T(j)-Te(j))^2/jmax)^(1/2) 
end 
 
function[Ue,xe,te]=truesolution(a,f,g1,g2,xspan,tspan,ngrid) 
n=ngrid(1);%t份数 
m=ngrid(2);%x份数 
Ue=zeros(ngrid); 
xe=linspace(xspan(1),xspan(2),m);%画网格  
te=linspace(tspan(1),tspan(2),n);%画网格 
for j=2:n 
    for i=2:m-1 
        for g=1:m-1 
Ue(j,i)=100-(400/(2*g-1)/pi)*sin((2*g-1)*pi*xe(j))*exp(-a*(2*g-1)^2*pi^2*te(i)) 
        end 
    end 
end 
 
function [U,x,t]=pdesolution(a,f,g1,g2,xspan,tspan,ngrid) 
n=ngrid(1);%t份数 
m=ngrid(2);%x份数 
h=range(xspan)/(m-1);%x网格长度  
x=linspace(xspan(1),xspan(2),m);%画网格  
k=range(tspan)/(n-1); %t网格长度 
t=linspace(tspan(1),tspan(2),n);%画网格 
U=zeros(ngrid); 
U(:,1)=g1(t);%边界条件 
U(:,m)=g2(t); 

  U(1,:)=f(x);%初值条件  
  %差分计算  
  for j=2:n  
      for i=2:m-1 
          U(j,i)=(1-2*a*k/h^2)*U(j-1,i)+a*k/h^2*U(j-1,i-1)+a*k/h^2*U(j-1,i+1);  
      end 
  end

