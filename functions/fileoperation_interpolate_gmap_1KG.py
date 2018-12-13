###Read in a 1KG genetic map in format:
#position rate(cM/Mb) genetic_map_position(cM)
###And a position file in format:
#pos1
#pos2
###And use linear interpolation to write out a new genetic map in formate
#pos1 cM
#pos2 cM

###GSJ 18/11/2017

import gzip
import numpy as np

def fileoperation_interpolate_gmap_1KG(mapfile, posfile, outfile, g_map_pos_IDx = 0, g_map_map_IDx = 2, g_map_delimiter = ' '):
    map_open = gzip.open if mapfile[-3:] == '.gz' else open
    gmap = []
    with map_open(mapfile, 'rb') as f_map:
        f_map.readline() #headers
        for line in f_map:
            split_line = line[0:-1].split(g_map_delimiter)
            #print g_map_pos_IDx, g_map_map_IDx
            #print split_line
            gmap.append([split_line[g_map_pos_IDx], split_line[g_map_map_IDx]])
    gmap = np.array(gmap, dtype = float)
    if gmap[0][0] == gmap[-1][0]:
        raise RuntimeError("Suspicious gmap behaviour - first position is same as last position, possible when IDXs are misspecified and chromosomes are encoded as integers")
    if len(np.unique(gmap[::,0])) <= 25:
        print "WARNING: There are <= 25 unique gmap positions. This could indicate that the position column is actually chromosome integers"
    if np.sum((gmap[::,1] - np.roll(gmap[::,1], 1)) < 0) != 1:
        raise RuntimeError("gmap cM coordinates appear not to be always increasing or equal as position increases?")
    #The positions and cM IDx are specified as g_map_pos_IDx and g_map_map_IDx
    pos_open = gzip.open if posfile[-3:] == '.gz' else open
    pos = []
    with pos_open(posfile, 'rb') as f_pos:
        for line in f_pos:
            pos.append(line[0:-1])
    pos = np.array(pos, dtype = int)
    #Intepolate.
    curr_gmap_pos = 0
    max_pos = len(gmap) - 1
    interp_map = []
    for curr_pos in pos:
        while gmap[curr_gmap_pos][0] < curr_pos:
            if curr_gmap_pos == max_pos:
                #Escape the loop if at end of the gmap.
                break
            else:
                curr_gmap_pos += 1
        if curr_gmap_pos == max_pos:
            #curr_pos is after the end of the gmap. Keep on adding this gmap position as no information.
            interp_map.append(gmap[curr_gmap_pos][1])
        elif curr_gmap_pos == 0:
            #curr_pos is before the start of the gmap. Keep on adding this gmap position as no information.
            interp_map.append(gmap[curr_gmap_pos][1])
        elif gmap[curr_gmap_pos][0] == curr_pos:
            interp_map.append(gmap[curr_gmap_pos][1])
        else:
            prev_curr_pos = [gmap[curr_gmap_pos - 1][0], gmap[curr_gmap_pos][0]]
            rel_dist = (curr_pos - prev_curr_pos[0]) / float(prev_curr_pos[1] - prev_curr_pos[0])
            interp_map.append(gmap[curr_gmap_pos - 1][1] + (rel_dist * (gmap[curr_gmap_pos][1] - gmap[curr_gmap_pos - 1][1])))
    #Write out
    with open(outfile, 'wb') as f_out:
        for i in range(len(pos)):
            f_out.write('\t'.join(['%d' %(pos[i]), '%.12f' %(interp_map[i])]) + '\n')
    return None
