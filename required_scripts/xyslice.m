function myslice = yzslice(f,kk)
% Returns the kkth yz slice of a 3D spins field
myslice=squeeze(f(kk,:, :));
