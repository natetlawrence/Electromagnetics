function [result Q]=Cyl2D_ML_PW_v12(Lambda,Radius,epsilon,isEz)
%v02 will work for all 3 dipole orientations
%   calculation of the E field for Hz polarization will no longer be done
%   using numerical derivative but instead an analytical expression for the
%   derivative is derived and scattered field coeffcients are used to
%   calculate field directly
%   
%   some code is included for Qscat and Qext but Im not sure its correct
%   for point source
%
%v2.2 will use direct calculation of source fields and Er Ephi from scat
%   coefficients, also includes Qscat and Qext 
%
%code will find the fields from a multilayered cylinder in 2D with a plane
%wave excitation at wavelength Lambda
%Radius is vector giving the thicknesses of each layer of the cylinder
%Epsilon gives the electric permittivity of each layer, must be 1 longer
%than Radius
%isEz gives the polarization of the problem, isEz=1, means E in z direciton
%result is Ez if isEz=1
%result is [Ex Ey] if isEz=2
clear global hnk jnk hpnk jpnk hnk1 jnk1 hpnk1 jpnk1 
global hnk jnk hpnk jpnk hnk1 jnk1 hpnk1 jpnk1 

epsilon=materials(Lambda,epsilon);
nParticle=sqrt(epsilon);

M=length(Radius);
k=2*pi/Lambda.*nParticle;

c=3e8;                      %speed light in m/s
omega=2*pi*c/(Lambda*1e-9);        %angular frequency
e0=8.85e-12;
imp_freespace=376.73;
E0=e0*1e-9*omega*imp_freespace; %constant required to normalize Er,Ephi when doing Hz polarization

%number of modes to be plotted
nmax=real(ceil(max(k)*max(Radius)+4*(max(k)*max(Radius))^(1/3)+2));
nmax=30;
%points in the where field will be plotted
plot_lims=250;

sz=3;
npts=250;
numx=npts;   %determines x values that field will be plotted at
xstart=-sz*max(Radius);
xstop=sz*max(Radius);
%xstart=-plot_lims;
%xstop=plot_lims;
x=xstart:(xstop-xstart)/(numx-1):xstop;

numy=npts;   %determines y values that field will be plotted at
ystart=-sz*max(Radius);
ystop=sz*max(Radius);
%ystart=-plot_lims;
%ystop=plot_lims;
y=ystart:(ystop-ystart)/(numy-1):ystop;

[X Y Z]=meshgrid(x,y,0);    %creates 2-D grid of values where each point
%is the x or y value at that point
[phi r z]=cart2pol(X,Y,Z);  %converts values to cylindrical coordinates

%matrix telling whether source is in each layer
%k_in will be the k vector of each point for field plotting
r_in=zeros(M,size(r,1),size(r,2));
k_in=zeros(size(r));
eps_in=zeros(size(r));
for kk=1:M
    r_in(kk,:,:)=r<Radius(kk);
end
r_in(M+1,:,:)=ones(size(r,1),size(r,2))-squeeze(r_in(M,:,:));
r_in(2:M,:,:)=r_in(2:M,:,:)-r_in(1:M-1,:,:);
for kk=1:M+1
    k_in=k_in+squeeze(r_in(kk,:,:)).*k(kk);
    eps_in=eps_in+squeeze(r_in(kk,:,:)).*epsilon(kk);
end

%creates matricies that are stacks of original matrix with height nmax
for nn=1:2*nmax+1
    phi_nmax(nn,:,:)=phi;
    r_nmax(nn,:,:)=r;
    n_nmax(nn,:,:)=(nn-nmax-1).*ones(size(r,1),size(r,2));
    eps_in_nmax(nn,:,:)=eps_in;
end
% phi_nmax=squeeze(phi_nmax);
% r_nmax=squeeze(r_nmax);
% n_nmax=squeeze(n_nmax);


for nn=-nmax-1:nmax+1
    hnk(:,nn+1+nmax+1)=besselh(nn,k(1:(length(k)-1)).*Radius);
    jnk(:,nn+1+nmax+1)=besselj(nn,k(1:(length(k)-1)).*Radius);
    hnk1(:,nn+1+nmax+1)=besselh(nn,k(2:(length(k))).*Radius);
    jnk1(:,nn+1+nmax+1)=besselj(nn,k(2:(length(k))).*Radius);
end
for nn=-nmax:nmax
    hpnk(:,nn+1+nmax)=transpose((k(1:(length(k)-1))./2)).*(hnk(:,nn+1+nmax)-hnk(:,nn+1+nmax+2));
    jpnk(:,nn+1+nmax)=transpose((k(1:(length(k)-1))./2)).*(jnk(:,nn+1+nmax)-jnk(:,nn+1+nmax+2));
    hpnk1(:,nn+1+nmax)=transpose((k(2:(length(k)))./2)).*(hnk1(:,nn+1+nmax)-hnk1(:,nn+1+nmax+2));
    jpnk1(:,nn+1+nmax)=transpose((k(2:(length(k)))./2)).*(jnk1(:,nn+1+nmax)-jnk1(:,nn+1+nmax+2));
end
hnk=hnk(:,2:size(hnk,2)-1);
jnk=jnk(:,2:size(jnk,2)-1);
hnk1=hnk1(:,2:size(hnk1,2)-1);
jnk1=jnk1(:,2:size(jnk1,2)-1);

if isEz
    SourceDir=3;
else 
    SourceDir=1;
end
ab=Cly2D_ML_coeff_v12(epsilon,Radius,nmax,Lambda,NaN,1,SourceDir);
Asn=zeros(size(r_nmax));
Bsn=zeros(size(r_nmax));
for nn=-nmax:nmax
    for mm=1:M+1
        Asn(nn+nmax+1,:,:)=squeeze(Asn(nn+nmax+1,:,:))+ab(mm,nn+nmax+1).*ones(size(r)).*squeeze(r_in(mm,:,:));
        Bsn(nn+nmax+1,:,:)=squeeze(Bsn(nn+nmax+1,:,:))+ab(mm+M+1,nn+nmax+1).*ones(size(r)).*squeeze(r_in(mm,:,:));
    end
end
Asn=squeeze(Asn);
Bsn=squeeze(Bsn);

%precalculated bessel functions for plotting fields
hnk2=zeros([2*nmax+1 size(r)]);
jnk2=zeros([2*nmax+1 size(r)]);
hpnk2=zeros([2*nmax+1 size(r)]);
jpnk2=zeros([2*nmax+1 size(r)]);
for nn=-nmax-1:nmax+1
    hnk2(nn+2+nmax,:,:)=besselh(nn,k_in.*r); %H_n(k_l*r) a set is calculated for each layer
    jnk2(nn+2+nmax,:,:)=besselj(nn,k_in.*r);
end
hpnk2=(permute(repmat(k_in/2,[1 1 2*nmax+1]),[3 1 2])).*(hnk2(1:size(hnk2,1)-2,:,:)-hnk2(3:size(hnk2,1),:,:));%H'_n(k_l*r) a set is calculated for each layer
jpnk2=(permute(repmat(k_in/2,[1 1 2*nmax+1]),[3 1 2])).*(jnk2(1:size(jnk2,1)-2,:,:)-jnk2(3:size(jnk2,1),:,:));

hnk2=hnk2(2:size(hnk2,1)-1,:,:);
jnk2=jnk2(2:size(jnk2,1)-1,:,:);
hnk2=squeeze(hnk2);
jnk2=squeeze(jnk2);

if isEz==1
    Ez_temp=(Asn.*jnk2+Bsn.*hnk2).*exp(-1i.*n_nmax.*phi_nmax);
    Ez_source=exp(1i.*k(length(k)).*X);
    Ez=squeeze(sum(Ez_temp,1));
    
    Ez=squeeze(Ez)+Ez_source.*squeeze(r_in(M+1,:,:));
    result(:,:,1)=Ez;
    result(:,:,3)=X;
    result(:,:,4)=Y;
elseif isEz==0
    E_r_temp=(n_nmax./(eps_in_nmax.*E0.*r_nmax)).*(Asn.*jnk2+Bsn.*hnk2).*exp(-1i.*n_nmax.*phi_nmax);
    E_phi_temp=(1./(1i*eps_in_nmax*E0)).*(Asn.*jpnk2+Bsn.*hpnk2).*exp(-1i.*n_nmax.*phi_nmax);

    E_r=squeeze(sum(E_r_temp,1)); %sums field over different modes n
    E_phi=squeeze(sum(E_phi_temp,1)); %sums field over different modes n
    
    Ex=cos(phi).*E_r-sin(phi).*E_phi;
    Ey=sin(phi).*E_r+cos(phi).*E_phi;

    Ey_source=exp(1i.*k(M+1).*X);
    
    Ey=Ey+Ey_source.*squeeze(r_in(M+1,:,:));
    result(:,:,1)=Ex;
    result(:,:,2)=Ey;
    result(:,:,3)=X;
    result(:,:,4)=Y;
end

%Calculates the scat and ext efficiencies from the field coefficients
%taken from BH p.204
%bn is the coefficients for the outgoing field of outermost layer
bn=squeeze(Bsn(:,1,1));
%this is old code for angular scattering but I dont know why!
for nn=-nmax:nmax
    bn(nn+nmax+1)=-bn(nn+nmax+1)./1i^nn;
end
Qext=bn(nmax+1);
for nn=1:nmax
    Qext=Qext+2*bn(nn+nmax+1);
end
Qext=real(Qext);
Qscat=abs(bn(nmax+1)).^2;
for nn=1:nmax
    Qscat=Qscat+2*abs(bn(nn+nmax+1))^2;
end
Qscat=Qscat.*2/(2.*pi./Lambda.*max(Radius));
Qext=Qext.*2./(2.*pi./Lambda.*max(Radius));
Q=[Qscat;Qext];

