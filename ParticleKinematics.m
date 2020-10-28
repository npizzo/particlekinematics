%% Kinematic behavior
cd /Users/npizzo/Documents/Research/ParticleKinematics/DLS_proj/
%
N = 512;
ML = 2*pi;
wl = 2*pi;
g = 1; 
xi = 0:1/N:1-1/N; 
x=xi;
% a=0.275;
a=1;
y=a*cos(2*pi*xi);
f=a*sin(2*pi*xi);
xi2= 0:max(xi)/(N-1):max(xi);
% xx=[x(1):X/(npts-1):x(end)];
xx=interp1(xi,x,xi2)';
clf
yy=interp1(x,y,xx);
ff=interp1(x,f,xx);
clf
plot(xx,yy)
hold on
plot(xx,ff)
x3=xx*2*pi;
% F=pt2(:,2);
%%
tl=15;
delete xc.txt
delete yc.txt 
delete fc.txt 
delete bw.txt
delete wl.txt
delete S.txt
% delete k.txt
BW=1;
% save xc.txt x_f -ascii
% save yc.txt y_f -ascii
% save fc.txt f_f -ascii
save xc.txt x3 -ascii
save yc.txt yy -ascii
save fc.txt ff -ascii
save S.txt a -ascii
save bw.txt BW -ascii
save tl.txt tl -ascii
save wl.txt wl -ascii 
% save k.txt k -ascii
unix('./run2.sh'); 
%%
formatSpec='%f'; 
fileID = fopen(['x.txt'], 'r');
xout = fscanf(fileID,formatSpec);
fclose(fileID); 
fileID = fopen(['y.txt'], 'r');
yout = fscanf( fileID, formatSpec);
fclose( fileID );
fileID = fopen(['phi.txt'], 'r');
pout = fscanf(fileID,formatSpec);
fclose(fileID);
fileID=fopen(['t.txt'], 'r');
to = fscanf(fileID,formatSpec);
fclose(fileID);
t=to(1:N:end);
nout_0=length(t);
nout=round(nout_0)-1;
data_o =zeros(N,nout_0);
x_o=zeros(N,nout_0);
y_o=zeros(N,nout_0);
p_o=zeros(N,nout_0);
for i=1:N
for j=1:nout_0
x_o(i,j)=xout(i+N*(j-1),1);
y_o(i,j)=yout(i+N*(j-1),1);
p_o(i,j)=pout(i+N*(j-1),1);
end
end
xi=0:2*pi/N:2*pi*(1-1/N);
dxi=abs(xi(2)-xi(1));
% reinterpolate
M=1*N;
xi2=0:2*pi/M:2*pi*(1-1/M);
x_1=zeros(M,length(t));
y_1=x_1;
for i=1:length(t)
    x_1(:,i)=interp1(xi,x_o(:,i),xi2);
    y_1(:,i)=interp1(xi,y_o(:,i),xi2);
end
t2=0:max(t)/1000:max(t);
x_2=zeros(M,length(t2));
y_2=x_2;
for i=1:M
    x_2(i,:)=interp1(t,x_1(i,:),t2);
    y_2(i,:)=interp1(t,y_1(i,:),t2);
end
%%
clf
set(gca,'fontsize',28)
for i=1:1:length(t)
plot(x_o(:,i),y_o(:,i),'.r')
xlim([0 2*pi])
ylim([-1 1]) 
pause(0.1)
hold off
end
% set(gca,'fontsize',28)
xlabel('X','interpreter','latex')
ylabel('$Y$','interpreter','latex')
%% compute velocity 
u=zeros(N,length(t));
v=zeros(N,length(t));
for i=1:N
   u(i,:)=gradient(x_o(i,:))./gradient(t');
   v(i,:)=gradient(y_o(i,:))./gradient(t');
end 
% compute accelerations 
ax=u; ay=u;
for i=1:N
    ax(i,:)=gradient(u(i,:))./gradient(t');
    ay(i,:)=gradient(v(i,:))./gradient(t');
end
%% find phase speed of max, max slope, max curvature
m = zeros(length(t),1);
ind = m;
xm = m; mx = m; indx = m;
xsm = xm; um = m; usm = m;
for i=1:length(t)
    
    [m(i),ind(i)]=max(y_o(:,i));
    
    xm(i,1)=x_o(ind(i),i);
    
    um(i,1)=u(ind(i),i);
    
    [mx(i), indx(i)] = min((gradient(y_o(:,i))./...
    gradient(x_o(:,i))));

    xsm(i,1)=x_o(indx(i),i);
    
    usm(i,1)=u(indx(i),i);
end
%% visualize this
clf
for i=1:length(t)
plot(x_o(:,i),y_o(:,i),'k')
hold on
plot(x_o(ind(i),i),y_o(ind(i),i),'or')
hold on
plot(x_o(indx(i),i),y_o(indx(i),i),'ob')
hold off
ylim([-1 1]) 
xlim([0 2*pi])
pause(0.1)
% hold off
end
%%
clf
plot(t,gradient(unwrap(xm)))
hold on
plot(t,gradient(xm),'r')

%% look at the ratios there
cm=gradient(smooth(unwrap(xm),5))./gradient(t);
cmx=gradient(smooth(unwrap(xsm),5))./gradient(t);
clf
p1=plot(t,smooth(um./cm,5),'r','linewidth',3);
hold on
p2=plot(t,smooth(usm./cmx,5),'k','linewidth',3);
ylim([-1 3])
xlim([0 4])
set(gca,'fontsize',22)
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
l1=legend([p1 p2],'$u(max(\eta))/c(max(\eta))$',...
    '$u(max(|\eta_x|))/c(max(|\eta_x|))$');
set(l1,'interpreter','latex')