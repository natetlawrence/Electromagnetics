function [result]=Cly2D_ML_coeff_v12(epsilon,r,nmax,Lambda,source,isPW,SourceDir)
%v02 will work for all 3 dipole orientations
%SourceDir is the direciton of the source, 1->x, 2->y, 3->z
global hnk jnk hpnk jpnk hnk1 jnk1 hpnk1 jpnk1 jnk_s hnk_s

if not(isPW)
[rsource(2) rsource(1)]=cart2pol(source(1),source(2)); %source location in polar
rs=rsource(1);
end
if SourceDir==3
    isEz=1;
else
    isEz=0;
end

M=length(r); %number of layers

index=sqrt(epsilon);
k=index*2*pi/Lambda;

%to change polarization just have to chage the boundary conditions being
%satisfied, if e_nmax=epsilon.^(3/2), TM (E-field in z-dir), if e_nmax=epsilon.^(-1/2),
%TE (H-field in z-dir)

if isEz==1
    pol=0;
else
    pol=-1;
end

for nn=1:2*nmax+1
    e_nmax(:,nn)=epsilon.^(pol);
end
e_nmax=squeeze(e_nmax);

%functions which will be used in calculation
if not(isnan(source))
rtemp=[r rs];
rtemp=sort(rtemp);
slayer=find(rtemp==rs);
end

%now construct the matrix equation which has to be solved
%examples fo the matrix can be found in 
%multilayered cylinder derivation.doc
unknowns=zeros(1,2*M,2*nmax+1);
knowns=zeros(1,2*M,2*nmax+1);
matrix=zeros(2*M,2*M,2*nmax+1);

for ii=1:M
    if ii==1
        matrix(2*ii-1,2*ii-1,:) = jnk(ii,:);
        matrix(2*ii-1,2*ii,:) = -jnk1(ii,:);
        matrix(2*ii-1,2*ii+1,:) = -hnk1(ii,:);
        matrix(2*ii,2*ii-1,:) = e_nmax(ii,:).*jpnk(ii,:);
        matrix(2*ii,2*ii,:) = -e_nmax(ii+1,:).*jpnk1(ii,:);
        matrix(2*ii,2*ii+1,:) = -e_nmax(ii+1,:).*hpnk1(ii,:);
    elseif ii==M
        matrix(2*ii-1,2*ii-2,:) = jnk(ii,:);
        matrix(2*ii-1,2*ii-1,:) = hnk(ii,:);
        matrix(2*ii-1,2*ii,:) = -hnk1(ii,:);
        matrix(2*ii,2*ii-2,:) = e_nmax(ii,:).*jpnk(ii,:);
        matrix(2*ii,2*ii-1,:) = e_nmax(ii,:).*hpnk(ii,:);
        matrix(2*ii,2*ii,:) = -e_nmax(ii+1,:).*hpnk1(ii,:);
    else
        matrix(2*ii-1,2*ii-2,:) = jnk(ii,:);
        matrix(2*ii-1,2*ii-1,:) = hnk(ii,:);
        matrix(2*ii-1,2*ii,:) = -jnk1(ii,:);
        matrix(2*ii-1,2*ii+1,:) = -hnk1(ii,:);
        matrix(2*ii,2*ii-2,:) = e_nmax(ii,:).*jpnk(ii,:);
        matrix(2*ii,2*ii-1,:) = e_nmax(ii,:).*hpnk(ii,:);
        matrix(2*ii,2*ii,:) = -e_nmax(ii+1,:).*jpnk1(ii,:);
        matrix(2*ii,2*ii+1,:) = -e_nmax(ii+1,:).*hpnk1(ii,:);
    end
end

if M==1
    matrix=zeros(2*M,2*M,2*nmax+1);
    matrix(1,1,:) = jnk;
    matrix(1,2,:) = -hnk1;
    matrix(2,1,:) = e_nmax(1,:).*jpnk;
    matrix(2,2,:) = -e_nmax(2,:).*hpnk1;
end

%make the correct vector of the known C quantities
%derivation of the Csn coefficients can be found in
%source coefficient derivation.doc
if not(isPW)
    n=(-nmax:nmax);
    if SourceDir==3 %z-oriented dipole
        C_in=(hnk_s.*(1/(4*1i))).*exp(-1i.*n.*rsource(2));
        C_out=(jnk_s.*(1/(4*1i))).*exp(-1i.*n.*rsource(2));
    elseif SourceDir==1 %x-oriented dipole
        C_out=((-1/(8*1i)).*sqrt(epsilon(slayer)).*(jnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))+jnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
        C_in=((-1/(8*1i)).*sqrt(epsilon(slayer)).*(hnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))+hnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
    elseif SourceDir==2 %y-oriented dipole
        C_out=((1/(8)).*sqrt(epsilon(slayer)).*(jnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))-jnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
        C_in=((1/(8)).*sqrt(epsilon(slayer)).*(hnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))-hnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
    end
    if slayer==1 %source in first layer
        knowns(1,1,:) = -(C_out).*hnk(slayer,:);
        knowns(1,2,:) = -(C_out).*e_nmax(slayer,:).*hpnk(slayer,:);
    elseif slayer==M+1 %source in last layer
        knowns(1,2*M-1,:) = (C_in).*jnk1(slayer-1,:);
        knowns(1,2*M,:) = (C_in).*e_nmax(slayer,:).*jpnk1(slayer-1,:);
    else
        knowns(1,2*slayer-3,:) = (C_in).*jnk1(slayer-1,:);
        knowns(1,2*slayer-2,:) = (C_in).*e_nmax(slayer,:).*jpnk1(slayer-1,:);
        knowns(1,2*slayer-1,:) = -(C_out).*hnk(slayer,:);
        knowns(1,2*slayer,:) = -(C_out).*e_nmax(slayer,:).*hpnk(slayer,:);
    end
else
    knowns=zeros(1,2*M,2*nmax+1);
    for nn=-nmax:nmax
        knowns(1,2*M-1,nn+nmax+1) = 1i^nn.*jnk1(M,nn+nmax+1);
        knowns(1,2*M,nn+nmax+1) = 1i^nn.*e_nmax(M+1,nn+nmax+1).*jpnk1(M,nn+nmax+1);
    end
end

%solve the system for the unknown coefficients
for nn=1:2*nmax+1
    unknowns(1,:,nn)=matrix(:,:,nn)\transpose(knowns(1,:,nn));
    %unknowns(1,:,nn)=pinv(matrix(:,:,nn))*transpose(knowns(1,:,nn));
end

Aln=zeros(M+1,2*nmax+1);
Bln=zeros(M+1,2*nmax+1);

Aln(1,:)=unknowns(1,1,:);
Bln(1,:)=zeros(size(unknowns(1,1,:)));
Aln(M+1,:)=zeros(size(unknowns(1,2*M,:)));
Bln(M+1,:)=unknowns(1,2*M,:);

for mm=2:M
    Aln(mm,:)=unknowns(1,2*mm-2,:);
    Bln(mm,:)=unknowns(1,2*mm-1,:);
end

result=[Aln;Bln];



