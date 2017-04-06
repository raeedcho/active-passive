function tesselateFR(fr,thv,thf)
% create patched voronoi diagram with given firing rate

if size(fr,2)>1
    error('Too many firing rates, only one at a time')
end

% caxis([min(fr) max(fr)]);

[v,c] = voronoin([thf thv]);

% fr_offset = min(fr);
% fr_scale = max(fr)-min(fr);
% norm_fr = (fr-fr_offset)./fr_scale;

for i = 1:length(c)
    if all(c{i}~=1)
        patch(v(c{i},1),v(c{i},2),fr(i))
    end
end

set(gca,'xlim',[-pi pi],'ylim',[-pi pi])
colorbar