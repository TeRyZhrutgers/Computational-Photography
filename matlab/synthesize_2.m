function [imout] = synthesize_2(imin, sizeout, tilesize, overlap)
% function [imout] = synthesize(imin, sizeout, tilesize, overlap)


imin = double(imin);

imout = zeros(sizeout);
sizein = size(imin);

threshold = 0.95;
scale = 0.5;

for i=1:tilesize-overlap:sizeout(1)-tilesize+1,
  for j=1:tilesize-overlap:sizeout(2)-tilesize+1,

    if (i > 1) && (j > 1),
    % extract top shared region
      shared = imout(i:i+overlap-1,j:j+tilesize-1);
      sim_top = xcorr2(imin, shared);
      
      % trim invalid data at edges, and off bottom where we don't want
      % tiles to go over the edge
      sim_top = sim_top(overlap:end-tilesize+1,tilesize:end-tilesize+1);

      % extract left shared region, skipping bit already checked
      shared = imout(i+overlap:i+tilesize-1,j:j+overlap-1);
      sim_left =  xcorr2(imin, shared);
      sim_left = sim_left(tilesize:end-tilesize+overlap+1, overlap:end-tilesize+1);
      % sum(shared(:).^2); trim invalid data at edges, and where we
      % don't want tiles to go over the edges
      sim = sim_top + sim_left-scale*(abs(sim_top-sim_left));

      [ibest, jbest] = find(sim >=threshold*max(sim(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      over_top = mincut(imout(i:i+overlap-1,j:j+tilesize-1),imin(pos(1):pos(1)+overlap-1,pos(2):pos(2)+tilesize-1));
      imout(i:i+overlap-1,j:j+tilesize-1) = over_top;
      
      over_left = mincut(imout(i:i+tilesize-1,j:j+overlap-1),imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+overlap-1));
      imout(i:i+tilesize-1,j:j+overlap-1) = over_left;
      
      imout(i+overlap:i+tilesize-1,j+overlap:j+tilesize-1) = imin(pos(1)+overlap:pos(1)+tilesize-1,pos(2)+overlap:pos(2)+tilesize-1);
    elseif i > 1
      shared = imout(i:i+overlap-1,j:j+tilesize-1);
      sim_top = xcorr2(imin, shared);
      
      % trim invalid data at edges, and off bottom where we don't want
      % tiles to go over the edge
      sim_top = sim_top(overlap:end-tilesize+1,tilesize:end-tilesize+1);

      [ibest, jbest] = find(sim_top >= threshold*max(sim_top(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      over = mincut(imout(i:i+overlap-1,j:j+tilesize-1),imin(pos(1):pos(1)+overlap-1,pos(2):pos(2)+tilesize-1));
      imout(i:i+overlap-1,j:j+tilesize-1) = over;
      
      imout(i+overlap:i+tilesize-1,j:j+tilesize-1) = imin(pos(1)+overlap:pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
    elseif j > 1
      shared = imout(i+overlap:i+tilesize-1,j:j+overlap-1);
      sim_left =  xcorr2(imin, shared);
      % sum(shared(:).^2); trim invalid data at edges, and where we
      % don't want tiles to go over the edges
      sim_left = sim_left(tilesize:end-tilesize+overlap+1, overlap:end-tilesize+1);

      [ibest, jbest] = find(sim_left >= threshold*max(sim_left(:)));
      c = ceil(rand * length(ibest));
      pos = [ibest(c) jbest(c)];
      over = mincut(imout(i:i+tilesize-1,j:j+overlap-1),imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+overlap-1));
      imout(i:i+tilesize-1,j:j+overlap-1) = over;
      imout(i:i+tilesize-1,j+overlap:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2)+overlap:pos(2)+tilesize-1);
    else
      pos = ceil(rand([1 2]) .* (sizein-tilesize+1));
      imout(i:i+tilesize-1,j:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
    end


   % imout(i:i+tilesize-1,j:j+tilesize-1) = imin(pos(1):pos(1)+tilesize-1,pos(2):pos(2)+tilesize-1);
  end
end

function [over] = mincut(old,new)

size_in = size(old);
if size_in(1) < size_in(2),
    old = old';
    new = new';
end
size_now = size(old);
distance = (old-new).^2;

row = size_now(1);
col = size_now(2);
for i=2:1:row,
    distance(i,1) = distance(i,1)+min(distance(i-1,1),distance(i-1,2));
    for j=2:1:col-1,
        distance(i,j) = distance(i,j)+min(distance(i-1,j-1),min(distance(i-1,j),distance(i-1,j+1)));
    end
    distance(i,col) = distance(i,col)+min(distance(i-1,col-1),distance(i-1,col));
end
  
over = zeros(size_now);


[data,index]  = min(distance(row,:));
over(row,1:index) = old(row,1:index);
if(index <col),
    over(row,index+1:col) = new(row,index+1:col);
end
for i=row-1:-1:1,
    if(index==1),
        [data,index] = min(distance(i,index:index+1));
    elseif(index==col),
        [data,index] = min(distance(i,index-1:index));
    else
        [data,index] = min(distance(i,index-1:index+1));
    end
    
    over(i,1:index) = old(i,1:index);
    if(index <col),
        over(i,index+1:col) = new(i,index+1:col);
    end
    
end

if size_in(1)<size_in(2),
    over = over';
end


