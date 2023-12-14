
clear all; tstart=tic; k=100; t=0.01; kt=k*t; qin=100;

ny=32;ny=min(ny,90);ny=round(ny/2)*2;nx=ny*2;nel=nx*ny;

nys=ny/2;Ntsa=nys*(nx-1)+nx*(nys+1) ;Ntca=(nys+1)*(nx-1);
Kca=kt/6*[4 -1 -2 -1;-1 4 -1 -2;-2 -1 4 -1;-1 -2 -1 4];

lc=zeros(nx*nys,4);my=nys+1;
for i=1:nx;
for j=1:nys;
lc((i-1)*nys+j,1)=my*(i-1)+j;
lc((i-1)*nys+j,2)=lc((i-1)*nys+j,1)+my;
lc((i-1)*nys+j,3)=lc((i-1)*nys+j,2)+1;
lc((i-1)*nys+j,4)=lc((i-1)*nys+j,1)+1;
end;
end;
K=zeros(Ntca+2*(nys+1),Ntca+2*(nys+1));
for n=1:nx*nys;
for i=1:4;
for j=1:4;
K(lc(n,i),lc(n,j))=K(lc(n,i),lc(n,j))+Kca(i,j);
end;
end;
end
Ksca=K(nys+2:(nx+1)*(nys+1)-nys-1,nys+2:(nx+1)*(nys+1)-nys-1); 
gca=zeros(Ntca,1);
for i=1+nys:my:nys*2*my-my;gca(i)=2*qin/nx;end;gca(ny*my)=qin/nx ;
tca=Ksca\gca;disca=gca'*tca;clear Ksca;clear Kca;

Ksa=kt/2*[5 1 -3 -3; 1 5 -3 -3; -3 -3 5 1; -3 -3 1 5] ;
ls=zeros(nx*nys,4);
for j=1:nx;
for i=1:nys;
ls((j-1)*nys+i,1)=(j-1)*(ny+1)+nys+i+1;
ls((j-1)*nys+i,2)=ls((j-1)*nys+i,1)-1;
ls((j-1)*nys+i,3)=j*(ny+1)+i;
ls((j-1)*nys+i,4)=ls((j-1)*nys+i,3)-ny-1;
end;
end
K=zeros(Ntsa+ny,Ntsa+ny);
for n=1:nx*nys;
for i=1:4;
for j=1:4;
K(ls(n,i),ls(n,j))=K(ls(n,i),ls(n,j))+Ksa(i,j) ;
end;
end;
end;
Kssa=K(nys+1:Ntsa+nys,nys+1:Ntsa+nys);clear K; 
gsa=zeros(Ntsa,1);
for i=1+nys:2*nys+1:ny*(ny+1)-nys;gsa(i)=qin/(ny);end;
tsa=Kssa\gsa;dissa=gsa'*tsa;clear Kssa;

figure('Position',[1 1 1600 800]); axes('fontsize',15)
xnod=0:2/nx:5;px=1/nx:2/nx:5-1/nx ;
test=[0 (tca(1:my:(nx-2)*my+1))' (0:2/nx:.5)*0 (tca((nx-1)*my:-my:my))'...
(0:2/nx:.5)*0];
td=[tsa(1:ny+1:(2*nys+1)*nx-ny)' (0:2/nx:.5-2/nx)*0 ...
tsa((ny+1)*nx-nys :-(ny+1) :nys+1 )' (0:2/nx:.5-2/nx)*0];
plot(xnod,test,'k');hold on;
plot(px,td,'.k','MarkerSize',12);hold on; grid on;
axis([-.02 5.02 -max(tca)/100 max(max(tca),max(td))/.98])
Error=(dissa-disca)*200/(dissa+disca);telapsed=toc(tstart);
title([num2str(nx),' x ',num2str(ny),', CA: ',num2str(disca,'%0.6g'),...
', SA: ',num2str(dissa,'%0.6g'),', Er: ',num2str(Error,'%0.2g'),...
' %, ','cpu: ',num2str(telapsed,'%0.3g')],'fontsize', 15)
legend('CA ','SA ','Location','NorthWest')

max(test)

