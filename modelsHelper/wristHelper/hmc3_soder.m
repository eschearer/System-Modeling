function[T,rms]=hmc3_soder(P1,P2)
% Estimates rigid body pose T from marker coordinates
%
% Input:
%	P1		(N x 3 matrix) coordinates of N markers in coordinate system 1
%	P2		(N x 3 matrix) coordinates of the same markers, now measured in coordinate system 2
%
% Output:
%	T		(4x4 matrix) homogeneous transformation matrix that transforms P1 coordinates
%			to P2 coordinates (P2 = T*P1)
%	rms		(scalar, optional) rms fit error of the rigid body model: P2 = T*P1 + error
%
% Method:
%	Soderkvist I. and Wedin P.A. (1993) Determining the movements of the skeleton using
%	well-configured markers.  Journal of Biomechanics, 26:1473-1477.    

	% error checking
	[nmarkers,ndim1]=size(P1);
	[nmarkers2,ndim2]=size(P2);
	if (ndim1 ~= 3) || (ndim2 ~= 3)
		error('hmc3_soder: matrix with marker coordinates must have 3 columns (x,y,z)');
	end
	if (nmarkers2 ~= nmarkers)
		error('hmc3+soder: matrices with marker coordinates must have same number of rows');
	end
		
	% construct matrices A and B, subtract the mean so there is only rotation
	m1=mean(P1);
	m2=mean(P2);
	for i=1:nmarkers
	  A(i,:)=P1(i,:)-m1;
	  B(i,:)=P2(i,:)-m2;
	end
	A = A';
	B = B';

	% The singular value decomposition to calculate rotation matrix R
	C=B*A';
	[P,T,Q]=svd(C);
	R=P*diag([1 1 det(P*Q')])*Q';

	% Calculate the translation vector from the centroid of all markers
	d=m2'-R*m1';

	% combine them into a 4x4 matrix
	T = [R d ; 0 0 0 1];

	% calculate RMS value of residuals
	sumsq = 0;
	for i=1:nmarkers
	  P2model = R*P1(i,:)' + d;
	  sumsq = sumsq + norm(P2model-P2(i,:)')^2;
	end
	rms = sqrt(sumsq/3/nmarkers);
end