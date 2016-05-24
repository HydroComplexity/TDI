if (lev==1 || lev==num_lev)     
  figure(1)
  subaxis(1,2,num_fig,'MT',0.1,'ML',0.1,'MR',0.1,'SpacingHoriz',0.11,'SpacingVert',0.13)
  imagesc(Puddle); hold on
  colormap('gray')
  xlabel('X','FontSize', 12);
  ylabel('Y','FontSize', 12);
  title(['Depression Level ',num2str(lev)],'FontSize', 12);

  figure(2)
  subaxis(1,2,num_fig,'MT',0.1,'ML',0.1,'MR',0.1,'SpacingHoriz',0.11,'SpacingVert',0.13)
  imagesc(Threshold); hold on
  colormap('gray')
  xlabel('X','FontSize', 12);
  ylabel('Y','FontSize', 12);
  title(['Depression Level ',num2str(lev)],'FontSize', 12);    

  num_fig = num_fig + 1;
end