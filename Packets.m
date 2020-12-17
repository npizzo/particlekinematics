%% Focusing packets
cd /Users/npizzo/Documents/GitHub/particlekinematics/
N = 512;
ML = 100;
xo = -30;
x = xo:ML/N:xo+(ML)*(1-1/N);
wl = ML;
g = 9.81; 
% g=1;
k = (2*pi)^2/g; 
% k=2*pi;
w = sqrt(g*k); cg=g/2/w; 
tstart = cputime;
for ii = 2
D = 0.75;
xb = 15 / D;
tb = xb / cg;
S = 0.1 + (ii-1)*0.1;
M = 32;
kn = zeros(M,1);
wn = zeros(M,1);
an = zeros(M,1);
for i = 1 : M
    kn(i) = k*(1+D*(i-M/2)/M);
    wn(i) = sqrt(g*kn(i));
    an(i) = S/kn(i)/M;
end

eta1 = zeros(length(x),1);
for i = 1 : M
    eta1(:,1) = eta1(:,1) + an(i)*cos(kn(i)*(x'-xb)+wn(i)*tb);
end
clf
plot(x,eta1)
filt1=1/2*(tanh(0.25*(x--15))-tanh(0.25*(x-10)));
eta2=eta1.*filt1';
hold on
plot(x,eta2)
%
I=sqrt(-1);
Fs=1/(abs(x(2)-x(1)));

[out,freq]=positiveFFT(eta2',Fs);
eta=zeros(1,length(x)); 
    for j=1:length(x)
        eta(1,j) = sum(real(2*out(2:end).*exp(I*(2*pi*freq(2:end)*...
            (x(j)+30)))));
    end
eps2=sum((abs(2*out(2:end)).*2.*pi.*freq(2:end)));
phi=zeros(1,length(x)); 
    for j=1:length(x)
        phi(1,j) = sum(imag(2*sqrt(g)*out(2:end)./...
            sqrt(2*pi*(freq(2:end))).*...
           exp(I*((2*pi*freq(2:end))*...
            (x(j)+30)-sqrt(g*2*pi*freq(2:end))*(0-0*tb)))));
    end
% BW=powerbw(eta,2*pi*Fs,[],10)/k 
BW=D;
%
clf
x_f=(x)'-xo;
% x_f2=(x)'-xo+30;
y_f=eta;
f_f=phi;
plot(x_f,y_f,'k')
hold on
plot(x_f,f_f,'r')
%
% delete k.txt
save xc.txt x_f -ascii
save yc.txt y_f -ascii
save fc.txt f_f -ascii
save S.txt eps2 -ascii
save S0.txt S -ascii
save bw.txt BW -ascii
save bw0.txt BW -ascii
save C.txt BW -ascii 
save wl.txt wl -ascii 
% save k.txt k -as
unix('./run2.sh'); 
%
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
% need to create the structure 
% ii =1;
data.x = x_o;
data.y = y_o;
data.p = p_o;
data.t = t; 
data.N = N;
data.S = S;
data.BW = D;
data.xb = xb;
data.tb = tb;
data.k = k;
data.omega = w;
data.wl = wl;
data.x_0 = x_f;
data.y_0 = y_f;
data.f_0 = f_f;
foc{ii,1} = data; 
save('foc.mat', 'foc')
ii
cputime-tstart
end 
% save(filename, sprintf('a%d',ii))
% save('foc' ,data ,datestr(now), '.mat');
% D='data';
% F=sprintf('%s_%s.mat',data, datestr(now));
% save(fullfile(D,F))
%%
clf
set(gca,'fontsize',28)
for i=1:10:length(t)
plot(x_o(:,i),y_o(:,i),'r')
xlim([0 2*pi])
ylim([-1e-2 1e-2]) 
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