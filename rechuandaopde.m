function rechuandaopde 
  %�����������ݣ�����t�ķ�Χ�Ҹ�����ĿҪ��ȡ����20000���������pdf�е��� 
  a=0.00045;%a��ȡֵ 
  xspan=[0 10];%x��ȡֵ��Χ 
  tspan=[0 5400];%t��ȡֵ��Χ 
  ngrid=[1000 100];%�ָ�ķ�����ǰ�����t��ģ��������x��� 
  f=@(x)37;%��ֵ 
  g1=@(t)75;%�߽�����һ 
  g2=@(t)75;%�߽������� 
  [T,x,t]=pdesolution(a,f,g1,g2,xspan,tspan,ngrid);%���������õĺ��� 
  [x,t]=meshgrid(x,t); 
  mesh(x,t,T);%��ͼ�����Ұ����������Ƹ�Ϊx��t��T 
  xlabel('x') 
  ylabel('t') 
  zlabel('T') 
  T%����¶Ⱦ��� 
  dt=tspan(2)/ngrid(1);%t���� 
  h3000=3000/dt; 

h9000=9000/dt; 
h15000=15000/dt;%3000,9000,15000�£��¶ȷֱ���T�������Щ�� 
T3000=T(h3000,:) 
T9000=T(h9000,:) 
T15000=T(h15000,:)%�������ʱ���µ��¶ȷֲ� 
%���ٶ�����ʱ���µ��¶�-�������߻�ͼ����ͼ�������άͼ�Ľ��� 
 
%�ȶ�������,����Ҷ������ 
dx=xspan(2)/ngrid(2);%x���� 
sta=4*a*dt/(dx^2)*(sin(pi/2))^2; 
if sta>0,sta<2 
    fprintf('\n%s\n','���ȶ���') 
else  
    fprintf('\n%s\n','û���ȶ���') 
    error 
end 
 
%��ʵֵ���� 
[xe,te,Te]=truesolution(a,f,g1,g2,xspan,tspan,ngrid); 
[xe,te]=meshgrid(xe,te); 
mesh(xe,te,Te);%��ͼ�����Ұ����������Ƹ�Ϊxe��te��Te 
xlabel('xe') 
ylabel('te') 
zlabel('Te') 
Te%����¶Ⱦ��� 
 
%������ 
jmax=1/dx+1;%������� 
[rms]=wuchajisuan(T,Te,jmax) 
rms%������ 

 
function [rms]=wuchajisuan(T,Te,jmax) 
for j=1:jmax 
    rms=((T(j)-Te(j))^2/jmax)^(1/2) 
end 
 
function[Ue,xe,te]=truesolution(a,f,g1,g2,xspan,tspan,ngrid) 
n=ngrid(1);%t���� 
m=ngrid(2);%x���� 
Ue=zeros(ngrid); 
xe=linspace(xspan(1),xspan(2),m);%������  
te=linspace(tspan(1),tspan(2),n);%������ 
for j=2:n 
    for i=2:m-1 
        for g=1:m-1 
Ue(j,i)=100-(400/(2*g-1)/pi)*sin((2*g-1)*pi*xe(j))*exp(-a*(2*g-1)^2*pi^2*te(i)) 
        end 
    end 
end 
 
function [U,x,t]=pdesolution(a,f,g1,g2,xspan,tspan,ngrid) 
n=ngrid(1);%t���� 
m=ngrid(2);%x���� 
h=range(xspan)/(m-1);%x���񳤶�  
x=linspace(xspan(1),xspan(2),m);%������  
k=range(tspan)/(n-1); %t���񳤶� 
t=linspace(tspan(1),tspan(2),n);%������ 
U=zeros(ngrid); 
U(:,1)=g1(t);%�߽����� 
U(:,m)=g2(t); 

  U(1,:)=f(x);%��ֵ����  
  %��ּ���  
  for j=2:n  
      for i=2:m-1 
          U(j,i)=(1-2*a*k/h^2)*U(j-1,i)+a*k/h^2*U(j-1,i-1)+a*k/h^2*U(j-1,i+1);  
      end 
  end

