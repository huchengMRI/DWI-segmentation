function seg_DWI(DWI_file, bval_file, mask_file, erode_mask)
%DWI_file: input DWI file name
%bval_file: input b-value file name
%mask_file: input mask of the DWI image
%erode_mask: a parameter for cleaning up voxels in the DWI image boundary
%because some "bad" voxels are usually included in nodif_brain_mask. The
%unit is voxel. e.g., 2
%in the example: seg_DWI('data.nii.gz', 'bvals', 'nodif_brain_mask.nii.gz', 2)

hdr = load_nii(DWI_file);
R = hdr.img;
bval = load(bval_file);
a = find(diff(sort(bval))>100);
[sbval, b] = sort(bval);
a1 = [1, a+1, length(bval)+1];


for i = 1:length(a)+1
    s{i} = b(a1(i):a1(i+1)-1);
end
[n1, n2, n3, n4] = size(R);
Sall = zeros(n1,n2,n3, length(a));
xb = zeros(1,length(a));
S0 = mean(squeeze(R(:,:,:,s{1})),4);
for i = 1:length(a)
    xb(i) = -log(mean(bval(s{i+1}))/mean(bval(s{2})));
    Sall(:,:,:,i) = mean(squeeze(R(:,:,:,s{i+1})),4);
end


alpha = zeros(size(S0));
beta = alpha;
h0 = load_nii(mask_file);
mask0 = h0.img;
[nx,ny,nz]= ind2sub([n1 n2 n3], find(mask0));
if (~exist('param.mat', 'file'))
for i = 1:length(nx)
    y = log(squeeze(Sall(nx(i),ny(i),nz(i), :))./S0(nx(i),ny(i),nz(i)));
    [a b] = polyfit(xb,y',1);
    alpha(nx(i),ny(i),nz(i)) = a(1);
    beta(nx(i),ny(i),nz(i)) = a(2);
end
save('param.mat', 'alpha', 'beta');
end
load param.mat;


nii.hdr = h0.hdr;    
nii.img = alpha.*double(h0.img)*1024;
save_nii(nii, 'alpha.nii');
nii.img = -beta.*double(h0.img)*1024;
save_nii(nii, 'beta.nii');
x = 0:0.001:4;
y = hist(-beta(beta<0),x);
[a,b] = max(y);
cutoff = 2*x(b);
csf = beta<-cutoff;
alpha_neg = 2-alpha;
%clean up voxels in the boundary, the level of cleanup nrv depends on nodif_brain_mask.
nrv = erode_mask;
[xx,yy,zz] = ndgrid(-nrv:nrv);
nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= nrv;
mask1 = imerode(mask0, nhood);
nii.img = alpha_neg.*(1-0.99*csf).*double(h0.img)*1024.*double(mask1);
nii.hdr = h0.hdr;
save_nii(nii, 'pseudo_T1w.nii');



