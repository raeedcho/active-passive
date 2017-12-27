function [ant_chans,post_chans] = split_array_Han(cds)
% splits array into anterior and posterior based on map file to check for
% area 2 vs area 5 responses.

chans = cat(1,cds.units.chan);
post_chans = unique(chans(cat(1,cds.units.rowNum) > 5));
ant_chans = unique(chans(cat(1,cds.units.rowNum) <= 5));