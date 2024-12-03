% T.Nikitopoulos
% Draw a molecular structure of Cyclic molecule
%
function drawmol(X)
 
global drawdata;
 
Nmol = length(X);
 
u = -X(1,:);   		% 1st atom always at (0,0,0)
for i=1:Nmol
   v = X(i,:);
   X(i,:) = v + u;   	% coords of other atoms while respecting
end                  	% distance constraints of the problem
 
X = [X ; X(1,:)];
 
handle=findobj('Type','figure','Name','DrawMol');
if isempty(handle);
 
  linewidth=4.6;   % original value 0.1^M
  pointwidth=15;  %    -"-    -"-   1^M
 
  figure('Name','DrawMol','NumberTitle','Off','BackingStore','Off','Color','k');
  set(gca,'DrawMode','Fast'); % alternative option of DrawMode used to be 'normal'
  set(gca,'position',[0.05 0.2 0.99 1.1],'units','normalized');
  set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w');
  axis('square');
  axis([-2, 2, -3, 3, -2, 2]);
  view(50,50);	% (45,15);
  hold on;
 
  for i=1:Nmol
     % evaluate X, Y, Z triple values for 3D
     Xpos = [X(i,1), X(i+1,1)];
     Ypos = [X(i,2), X(i+1,2)];
     Zpos = [X(i,3), X(i+1,3)];
 
     drawdata(i) = line(Xpos,Ypos,Zpos,'Color','r','LineStyle','-','LineWidth',linewidth,'EraseMode','xor');
  end % for i
 
  linewidth=1.0;   % original value 0.1
  Xpos = [X(2,1), X(8,1)];
  Ypos = [X(2,2), X(8,2)];
  Zpos = [X(2,3), X(8,3)];
 
  drawdata(Nmol+1) = line(Xpos,Ypos,Zpos,'Color','b','LineStyle','-','LineWidth',linewidth,'EraseMode','xor');
  drawdata(Nmol+2) = plot3(X(1,1),X(1,2),X(1,3),'Color','r','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
  drawdata(Nmol+3) = plot3(X(2,1),X(2,2),X(2,3),'Color','g','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
  drawdata(Nmol+4) = plot3(X(3,1),X(3,2),X(3,3),'Color','b','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
  drawdata(Nmol+5) = plot3(X(4,1),X(4,2),X(4,3),'Color','y','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
  drawdata(Nmol+6) = plot3(X(5,1),X(5,2),X(5,3),'Color','g','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
  drawdata(Nmol+7) = plot3(X(6,1),X(6,2),X(6,3),'Color','b','Marker','o','LineWidth',pointwidth,'EraseMode','xor');
 
  drawnow;
  return;
 
end % if isempty(handle);
 
for i=1:Nmol
   % evaluate X, Y, Z triple values for 3D^M
   Xpos = [X(i,1), X(i+1,1)];
   Ypos = [X(i,2), X(i+1,2)];
   Zpos = [X(i,3), X(i+1,3)];
 
   set(drawdata(i),'Xdata',Xpos,'Ydata',Ypos,'Zdata',Zpos);
   set(drawdata(i+1+Nmol),'Xdata',X(i,1),'Ydata',X(i,2),'Zdata',X(i,3));
end
 
Xpos = [X(2,1), X(8,1)];
Ypos = [X(2,2), X(8,2)];
Zpos = [X(2,3), X(8,3)];
set(drawdata(Nmol+1),'Xdata',Xpos,'Ydata',Ypos,'Zdata',Zpos);
 
figure(handle);
drawnow;
%
% end of drawmol
