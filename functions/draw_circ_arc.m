function [fig] = draw_circ_arc(center,radius,point1,point2,NameFigure,filled,num)

if abs(norm(center-point1)- norm(center-point2)) > 1e-14
    error("the distances between the points and the centers must be equal")
end

if nargin == 4
    NameFigure = 1;
    filled = 1;
    num = [2 20];
elseif nargin == 5
    filled = 1;
    num =[2 20];
elseif nargin == 6
    num = [2 20];
end
line1 = [linspace(center(1),point1(1),num(1))',linspace(center(2),point1(2),num(1))'];
line2 = [linspace(point2(1),center(1),num(1))',linspace(point2(2),center(2),num(1))'];

middle = (point2+point1)./2;
new_center = middle + radius.*(center-middle);

new_radius = norm(new_center-point1); % = |new_center-point2|
truc1 = point1 - new_center;
truc2 = point2 - new_center;

angl1 = atan(truc1(2)/truc1(1));
angl2 = atan(truc2(2)/truc2(1));

if truc1(1) < 0
    angl1 = angl1 + pi;
end
if truc2(1) < 0
    angl2 = angl2 + pi;
end

angl = linspace(angl1,angl2,num(2))';

arc = [new_radius.*cos(angl)+new_center(1),new_radius.*sin(angl)+new_center(2)];

fig = zeros(2*num(1)+num(2),2);
fig(1:num(1),:) = line1;
fig(num(1)+1:num(1)+num(2),:) = arc;
fig(num(1)+num(2)+1:2*num(1)+num(2),:) = line2;
if filled == 0
    figure(NameFigure)
    plot(fig(:,1),fig(:,2))
elseif filled == 1
    figure(NameFigure)
    fill(fig(:,1),fig(:,2),'blue','FaceAlpha',0.1)
end
end

