function tangle_visualize( out )
% input is the second output of tangleAnalysis.m

t1 = randsample(out.times_t1,1); % choose a random t1
thisT1 = out.allCombos == t1;

pairs2use = any(thisT1,2);
thisPair = out.allCombos(pairs2use,:);
q_T1 = out.q_allPairs(pairs2use,:);
X = out.X;
X_dot = out.X_dot;

%% Visualize the projections and how tangled they are

figure; set(gcf,'Color',[1 1 1]);
j = [];

% to make things nice to plot, set variance = 1
X = X/sqrt(sum(var(X)));
X_dot = X_dot/sqrt(sum(var(X_dot)));

step = 5;
for tt = 1:step:sum(pairs2use)
   
   
   delete(j)
   cla
   hold on;
   t2 = thisPair(tt,:) ~= t1;
   t2 = thisPair(tt,t2);
   
   % plot trajectories associated with the pair of times
   text(0, 1.5,'Trajectory through state space ("X")','FontSize',16,'HorizontalAl','center')
   text(0, -2.5,'Derivative of trajectory ("X_{dot}")','FontSize',16,'HorizontalAl','center')
   yOff = -4;
   plot(0,yOff,'rx','LineWidth',3)
   
   cond1 = out.conditionMask(t1); % which condition does t1 come from?
   cond2 = out.conditionMask(t2); % which condition does t2 come from?
   
   color1 = [.8 .2 .8];
   % color1 = 'b';
   plotProj(t1, 1, color1)
   if cond1 ~= cond2
      color2 = [.2 .8 .2];
%       color2 = 'k';
   else
      color2 = color1;
   end
   plotProj(t2, 2, color2)
   
   
   set(gca,'Visible','off')
   
   hightangling = 5e3;
   medtangling = 3e3;
   
   if q_T1(tt) > hightangling
      plotLine('r','r');
      
   elseif q_T1(tt) > medtangling
      plotLine([1 .5 0],[1 .5 0]);
      
   else % low tangling
      plotLine('y','k');
   end
   
   
   text(-2.2, yOff/2 +1, sprintf('t1: %1.0f, condition: %1.0f',t1, cond1),'FontSize',14,'HorizontalAl','left','Color','b')
   text(-2.2, yOff/2 +.8, sprintf('t2: %1.0f, condition: %1.0f',t2, cond2),'FontSize',14,'HorizontalAl','left')
   
   axis([-3 3 -6 2])
   pause
end


   function plotProj(timeIdx, timeID, projColor)
      cond = out.conditionMask(timeIdx);
      
      plot(0,0,'rx','LineWidth',3)
      plot(X(out.conditionMask == cond,1), X(out.conditionMask == cond,2),'Color',projColor)
      
      % plot derivatives (as arrows) associated with the pair of times
      if timeID == 1
         arrowColor = 'b';
         deriveMarker = 'bo';
         deriveLine = 'b-';
      else
         arrowColor = 'k';
         deriveMarker = 'ko';
         deriveLine = 'k-';
      end
      arrowsize = 8*norm(X_dot(timeIdx,[1:2]));
      try
         plot(X(timeIdx-10:timeIdx,1), X(timeIdx-10:timeIdx,2),'k','LineWidth',3)
         arrowMMC(X(timeIdx-1,1:2), X(timeIdx,1:2), [], arrowsize,[get(gca,'XLim') get(gca,'YLim')],arrowColor,'k');
      end
      
      % plot derivative trajectories
      plot(X_dot(out.conditionMask == cond,1), yOff+X_dot(out.conditionMask == cond,2),'Color',projColor)
      
      plot(X_dot(timeIdx,1),yOff+X_dot(timeIdx,2),deriveMarker,'MarkerFaceColor',arrowColor,'MarkerSize',10)
      plot([0 X_dot(timeIdx,1)],[yOff yOff+X_dot(timeIdx,2)],deriveLine,'LineWidth',2)
      
      
   end
%
   function plotLine(color, textColor)
      text(-2.5, yOff/2,sprintf('q(t1,t2): %1.2f',q_T1(tt)),'FontSize',14,'HorizontalAl','left','Color',textColor)
      plot([X_dot(t1,1) X_dot(t2,1)], yOff+[X_dot(t1,2) X_dot(t2,2)],'Color',color,'LineWidth',3)
      plot([X(t1,1) X(t2,1)], [X(t1,2) X(t2,2)],'Color',color,'LineWidth',3);
   end



end

