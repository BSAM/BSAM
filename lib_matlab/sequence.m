function[] = sequence(s0,nn1,nn2,step,var,varc,thres,toggle1,toggle2)
% s0 string with the name of the problem
% nn1 start frame
% nn2 end frame
% step of the representation
% All others as in 'interface3dcolor'

for nn=nn1:step:nn2
disp(sprintf('processing frame = %d \n', nn))
interface3dcolor(s0,nn,var,varc,thres,toggle1,toggle2)
end
