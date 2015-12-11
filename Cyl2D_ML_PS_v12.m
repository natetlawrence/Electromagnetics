function [result Gamma]=Cyl2D_ML_PS_v12(Lambda,Radius,epsilon,source,SourceDir)
%v02 will work for all 3 dipole orientations
%   calculation of the E field for Hz polarization will no longer be done
%   using numerical derivative but instead an analytical expression for the
%   derivative is derived and scattered field coeffcients are used to
%   calculate field directly
%   
%   some code is included for Qscat and Qext but Im not sure its correct
%   for point source
%v2.1 uses direct calculation of source fields instead of sum
%v2.2 will also calculate scattering cross section
%v2.3 calculates radiated power in a different way
%v2.4 calculate Gamma, Gamma=[G Gr Gnr]
%       G comes from LDOS, Gnr by integrating field, Gr from poynting
%       vector
%v2.5 corrects sign error of in field calculation which made Ex come out negative 

%code will find the fields from a multilayered cylinder in 2D with a point
%source excitation at wavelength Lambda
%Radius is vector giving the thicknesses of each layer of the cylinder
%Epsilon gives the electric permittivity of each layer, must be 1 longer
%than Radius
%result is Ez if isEz=1
%result is [Ex Ey] if isEz=2
clear global hnk jnk hpnk jpnk hnk1 jnk1 hpnk1 jpnk1 jnk_s hnk_s
global hnk jnk hpnk jpnk hnk1 jnk1 hpnk1 jpnk1 jnk_s hnk_s

[rsource(2) rsource(1)]=cart2pol(source(1),source(2)); %source location in polar
rs=rsource(1);

if SourceDir==3
    isEz=1;
else
    isEz=0;
end

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
nmax=100;
%points in the where field will be plotted
plot_lims=250;

sz=2.5;
npts=100;
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

%find layer which contains the source
rtemp=[Radius rs];
rtemp=sort(rtemp);
slayer=find(rtemp==rs);

%precalculated bessel funcitons for coefficient calculations
for nn=-nmax-1:nmax+1
    hnk(:,nn+1+nmax+1)=besselh(nn,k(1:(length(k)-1)).*Radius);
    jnk(:,nn+1+nmax+1)=besselj(nn,k(1:(length(k)-1)).*Radius);
    hnk1(:,nn+1+nmax+1)=besselh(nn,k(2:(length(k))).*Radius);
    jnk1(:,nn+1+nmax+1)=besselj(nn,k(2:(length(k))).*Radius);
    jnk_s(nn+1+nmax+1) = besselj(nn,k(slayer).*rs);
    hnk_s(nn+1+nmax+1) = besselh(nn,k(slayer).*rs);
end
for nn=-nmax:nmax
    hpnk(:,nn+1+nmax)=transpose((k(1:(length(k)-1))./2)).*(hnk(:,nn+1+nmax)-hnk(:,nn+1+nmax+2));
    jpnk(:,nn+1+nmax)=transpose((k(1:(length(k)-1))./2)).*(jnk(:,nn+1+nmax)-jnk(:,nn+1+nmax+2));
    hpnk1(:,nn+1+nmax)=transpose((k(2:(length(k)))./2)).*(hnk1(:,nn+1+nmax)-hnk1(:,nn+1+nmax+2));
    jpnk1(:,nn+1+nmax)=transpose((k(2:(length(k)))./2)).*(jnk1(:,nn+1+nmax)-jnk1(:,nn+1+nmax+2));
    hpnk_s(nn+1+nmax)=transpose((k(slayer)./2)).*(hnk_s(nn+1+nmax)-hnk_s(nn+1+nmax+2));
    jpnk_s(nn+1+nmax)=transpose((k(slayer)./2)).*(jnk_s(nn+1+nmax)-jnk_s(nn+1+nmax+2));
end
if SourceDir==3
    hnk_s=hnk_s(2:size(hnk_s,2)-1);
    jnk_s=jnk_s(2:size(jnk_s,2)-1);
end
hnk=hnk(:,2:size(hnk,2)-1);
jnk=jnk(:,2:size(jnk,2)-1);
hnk1=hnk1(:,2:size(hnk1,2)-1);
jnk1=jnk1(:,2:size(jnk1,2)-1);

ab=Cly2D_ML_coeff_v12(epsilon,Radius,nmax,Lambda,source,0,SourceDir);
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

Xs = X - ones(size(X))*source(1); %distance for each plot point from source
Ys = Y - ones(size(Y))*source(2);
[PHIs Rs] = cart2pol(Xs,Ys);

%sum field over multipole orders
if SourceDir==3
    %calculate fields form equation (4)
%contribution of scattered fields
E_z_temp=(Asn.*jnk2+Bsn.*hnk2).*exp(1i.*n_nmax.*phi_nmax);

Ez_source=besselh(0,k(slayer).*Rs)/(4*1i);

E_z=squeeze(sum(E_z_temp,1)); %sums field over different modes n

E_z=squeeze(E_z)+Ez_source.*squeeze(r_in(slayer,:,:));
result(:,:,1)=E_z;
result(:,:,3)=X;
result(:,:,4)=Y;

elseif SourceDir==1
%calculate fields form equation (4)
%contribution of scattered fields
%Er Ephi found in Hz pol derivation of e field
E_r_temp=-(n_nmax./(eps_in_nmax.*E0.*r_nmax)).*(Asn.*jnk2+Bsn.*hnk2).*exp(1i.*n_nmax.*phi_nmax);
E_phi_temp=(1./(1i*eps_in_nmax*E0)).*(Asn.*jpnk2+Bsn.*hpnk2).*exp(1i.*n_nmax.*phi_nmax);

Er_source=-1./(1i*E0*Rs).*(cos(PHIs)).*(-besselh(1,k(slayer)*Rs)/4);
Ephi_source=-k(slayer)/(2*1i*E0).*(sin(PHIs)/4).*(besselh(0,k(slayer)*Rs)-besselh(2,k(slayer)*Rs));

E_r=squeeze(sum(E_r_temp,1)); %sums field over different modes n
E_phi=squeeze(sum(E_phi_temp,1)); %sums field over different modes n

E_r=squeeze(E_r);
E_phi=squeeze(E_phi);
Ex=cos(phi).*E_r-sin(phi).*E_phi;
Ey=sin(phi).*E_r+cos(phi).*E_phi;

E_r_s=Er_source.*squeeze(r_in(slayer,:,:));
E_phi_s=Ephi_source.*squeeze(r_in(slayer,:,:));
Ex_s=cos(PHIs).*E_r_s-sin(PHIs).*E_phi_s;
Ey_s=sin(PHIs).*E_r_s+cos(PHIs).*E_phi_s;

result(:,:,1)=Ex+Ex_s;
result(:,:,2)=Ey+Ey_s;
result(:,:,3)=X;
result(:,:,4)=Y;

elseif SourceDir==2
%calculate fields form equation (4)
%contribution of scattered fields
%Er Ephi found in Hz pol derivation of e field
E_r_temp=-(n_nmax./(eps_in_nmax.*E0.*r_nmax)).*(Asn.*jnk2+Bsn.*hnk2).*exp(1i.*n_nmax.*phi_nmax);
E_phi_temp=(1./(1i*eps_in_nmax*E0)).*(Asn.*jpnk2+Bsn.*hpnk2).*exp(1i.*n_nmax.*phi_nmax);

Er_source=1./(1i*E0*Rs).*(-sin(PHIs)).*(-besselh(1,k(slayer)*Rs)/4);
Ephi_source=k(slayer)/(2*1i*E0).*(cos(PHIs)/4).*(besselh(0,k(slayer)*Rs)-besselh(2,k(slayer)*Rs));

E_r=squeeze(sum(E_r_temp,1)); %sums field over different modes n
E_phi=squeeze(sum(E_phi_temp,1)); %sums field over different modes n

E_r=squeeze(E_r);
E_phi=squeeze(E_phi);
Ex=cos(phi).*E_r-sin(phi).*E_phi;
Ey=sin(phi).*E_r+cos(phi).*E_phi;

E_r_s=Er_source.*squeeze(r_in(slayer,:,:));
E_phi_s=Ephi_source.*squeeze(r_in(slayer,:,:));
Ex_s=cos(PHIs).*E_r_s-sin(PHIs).*E_phi_s;
Ey_s=sin(PHIs).*E_r_s+cos(PHIs).*E_phi_s;

result(:,:,1)=Ex+Ex_s;
result(:,:,2)=Ey+Ey_s;
result(:,:,3)=X;
result(:,:,4)=Y;

end


%calculates radiated power Gr by integrating the opynting vector flux through
%a surface outside the structure. Same thing is done in BH, this can be
%calculated from the outgoing field coefficients
bn=squeeze(Bsn(:,1,1));
n=(-nmax:nmax);
if slayer==M+1  %if source is outside particle
    if SourceDir==1
        C_out=((-1/(8*1i)).*sqrt(epsilon(slayer)).*(jnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))+jnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
    elseif SourceDir==2
        C_out=((1/(8)).*sqrt(epsilon(slayer)).*(jnk_s(1:2*nmax+1).*exp(-1i.*(n-1).*rsource(2))-jnk_s(3:2*nmax+3).*exp(-1i.*(n+1).*rsource(2))));
    elseif SourceDir==3
        C_out=(jnk_s.*(1/(4*1i))).*exp(-1i.*n.*rsource(2));
    end
else
    C_out=zeros(size(bn));
end

Gr=0;
for nn=-nmax:nmax
    Gr=Gr+(bn(nn+nmax+1)+C_out(nn+nmax+1)).*conj(bn(nn+nmax+1)+C_out(nn+nmax+1));
end

%calculate the absorbed power in the material by numerically integrating
%the field
%img_eps give the imaginary part of epsilon for each poitn in the array
img_eps=imag(squeeze(eps_in));
for mm=1:M+1
    temp1=abs(result(:,:,1)).^2.*img_eps.*squeeze(r_in(mm,:,:));
    temp2=abs(result(:,:,2)).^2.*img_eps.*squeeze(r_in(mm,:,:));
    dArea=(x(2)-x(1))*(y(2)-y(1));
    Gnr(mm)=omega*e0*1e-9*sum(sum(temp1+temp2))*dArea;
end
Gnr=sum(Gnr);

%get rid of the extra calculated orders needed for derivatives
if not(SourceDir==3)
hnk_s=hnk_s(2:size(hnk_s,2)-1);
jnk_s=jnk_s(2:size(jnk_s,2)-1);
end

%calculate Gamma (G) from the LDOS
if SourceDir==3
    %A,B coefficients for source location
    An_s=ab(slayer,:);
    Bn_s=ab(slayer+M+1,:);
    
    E_z_temp=(An_s.*jnk_s+Bn_s.*hnk_s).*exp(1i.*n.*rsource(2));
    E_z=squeeze(sum(E_z_temp)); %sums field over different modes n

    %take imaginary part of greens funciton and add source contribution
    LDOS=imag(E_z+1/(4*1i));
    
elseif SourceDir==1
    %A,B coefficients for source location
    An_s=ab(slayer,:);
    Bn_s=ab(slayer+M+1,:);
    
    E_r_temp=-(n./(epsilon(slayer).*E0.*rsource(1))).*(An_s.*jnk_s+Bn_s.*hnk_s).*exp(1i.*n.*rsource(2));
    E_phi_temp=(1./(1i*epsilon(slayer)*E0)).*(An_s.*jpnk_s+Bn_s.*hpnk_s).*exp(1i.*n.*rsource(2));
    
    E_r=squeeze(sum(E_r_temp));
    E_phi=squeeze(sum(E_phi_temp));
    
    E_x=cos(rsource(2)).*E_r-sin(rsource(2)).*E_phi;
    E_y=sin(rsource(2)).*E_r+cos(rsource(2)).*E_phi;
    
    LDOS=imag(E_x+1/(8*1i));

elseif SourceDir==2
    %A,B coefficients for source location
    An_s=ab(slayer,:);
    Bn_s=ab(slayer+M+1,:);
    
    E_r_temp=-(n./(epsilon(slayer).*E0.*rsource(1))).*(An_s.*jnk_s+Bn_s.*hnk_s).*exp(1i.*n.*rsource(2));
    E_phi_temp=(1./(1i*epsilon(slayer)*E0)).*(An_s.*jpnk_s+Bn_s.*hpnk_s).*exp(1i.*n.*rsource(2));
    
    E_r=squeeze(sum(E_r_temp));
    E_phi=squeeze(sum(E_phi_temp));
    
    E_x=cos(rsource(2)).*E_r-sin(rsource(2)).*E_phi;
    E_y=sin(rsource(2)).*E_r+cos(rsource(2)).*E_phi;
    
    LDOS=imag(E_y+1/(8*1i));
end
G=LDOS;
Gamma=[G Gr Gnr];



