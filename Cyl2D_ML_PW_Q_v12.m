function [Q]=Cyl2D_ML_PW_Q_v12(Lambda,Radius,epsilon,isEz)
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

for nn=0:nmax+1
    hnk(:,nn+1+nmax+1)=besselh(nn,k(1:(length(k)-1)).*Radius);
    jnk(:,nn+1+nmax+1)=besselj(nn,k(1:(length(k)-1)).*Radius);
    hnk1(:,nn+1+nmax+1)=besselh(nn,k(2:(length(k))).*Radius);
    jnk1(:,nn+1+nmax+1)=besselj(nn,k(2:(length(k))).*Radius);
end
jnk(:,1:nmax+1)=fliplr(repmat((-1).^(1:(nmax+1)),[M 1]).*jnk(:,nmax+3:2*nmax+3));
hnk(:,1:nmax+1)=fliplr(repmat((-1).^(1:(nmax+1)),[M 1]).*hnk(:,nmax+3:2*nmax+3));
jnk1(:,1:nmax+1)=fliplr(repmat((-1).^(1:(nmax+1)),[M 1]).*jnk1(:,nmax+3:2*nmax+3));
hnk1(:,1:nmax+1)=fliplr(repmat((-1).^(1:(nmax+1)),[M 1]).*hnk1(:,nmax+3:2*nmax+3));
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


%Calculates the scat and ext efficiencies from the field coefficients
%taken from BH p.204
%bn is the coefficients for the outgoing field of outermost layer
bn=ab((M+1)*2,:);
%this is old code for angular scattering but I dont know why!
for nn=-nmax:nmax
    bn(nn+nmax+1)=-bn(nn+nmax+1)./1i^nn;
end
Qext=bn(nmax+1);
for nn=1:nmax
    Qext=Qext+2*bn(nn+nmax+1);
end
Qext=real(Qext);
Qscat=0;
for nn=-nmax:nmax
    Qscat=Qscat+abs(bn(nn+nmax+1))^2;
end
Qscat=Qscat.*2/(2.*pi./Lambda.*max(Radius));
Qext=Qext.*2./(2.*pi./Lambda.*max(Radius));
Q=[Qscat;Qext];

