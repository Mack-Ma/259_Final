% Rotating Descartes coordinates to a certain degree in a 2-D space
% Input raw data/the origin/rotation angle(-pi,pi,anticlockwise)
% Raw data is a n*2 matrix, the first colunm denotes x

function A=CRT_2D(data,origin,angle)
c=cosd(angle);
s=sind(angle);
x=data(:,1)-origin(1);
y=data(:,2)-origin(2);
x1=x*c+y*s;
y1=y*c-x*s;
A(:,1)=x1+origin(1);
A(:,2)=y1+origin(2);
end

% Memory Attention & Cognition(MAC) Lab %
% Edited by Tianye Ma 2018.3.31 %
