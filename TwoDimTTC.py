########################################################################################################
#
# Use function TTC(samples, 'dataframe') or TTC(samples, 'values') to compute two-dimensional Time-To-Collision.
#
#
# The first input is a pandas dataframe of vehicle pair samples, which should include the following columns.
# -----------------------------------------------------------------------------------
# x_i      :  x coordinate of the ego vehicle (usually assumed to be centroid)      |
# y_i      :  y coordinate of the ego vehicle (usually assumed to be centroid)      |
# vx_i     :  x coordinate of the velocity of the ego vehicle                       |
# vy_i     :  y coordinate of the velocity of the ego vehicle                       |
# hx_i     :  x coordinate of the heading direction of the ego vehicle              |
# hy_i     :  y coordinate of the heading direction of the ego vehicle              |
# length_i :  length of the ego vehicle                                             |
# width_i  :  width of the ego vehicle                                              |
# x_j      :  x coordinate of another vehicle (usually assumed to be centroid)      |
# y_j      :  y coordinate of another vehicle (usually assumed to be centroid)      |
# vx_j     :  x coordinate of the velocity of another vehicle                       |
# vy_j     :  y coordinate of the velocity of another vehicle                       |
# hx_j     :  x coordinate of the heading direction of another vehicle              |
# hy_j     :  y coordinate of the heading direction of another vehicle              |
# length_j :  length of another vehicle                                             |
# width_j  :  width of another vehicle                                              |
#------------------------------------------------------------------------------------
# The second input allows outputing a dataframe with inputed samples plus a new column named 'TTC', or mere TTC values.
# 
#
# If TTC==np.inf, the ego vehicle and another vehicle will never collide if they keep current speed.
# A negative TTC means the boxes of the ego vehicle and another vehicle are overlapping.
# This is due to approximating the space occupied by a vehicle with a rectangular.
# In other words, TTC<0 in this computation means the collision between the two vehicles almost (or although seldom, already) occurred.
#
########################## Copyright (c) 2022 Yiru Jiao <y.jiao-1@tudelft.nl> ###########################

# Import
import numpy as np
import warnings

# Functions
def line(point0, point1):
    x0, y0 = point0
    x1, y1 = point1
    a = y0 - y1
    b = x1 - x0
    c = x0*y1 - x1*y0
    return a, b, c

def intersect(line0, line1):
    a0, b0, c0 = line0
    a1, b1, c1 = line1
    D = a0*b1 - a1*b0 # D==0 then two lines overlap
    x = (b0*c1 - b1*c0)/D
    y = (a1*c0 - a0*c1)/D
    return np.array([x, y])

def ison(line_start, line_end, point):
    crossproduct = (point[1]-line_start[1])*(line_end[0]-line_start[0]) - (point[0]-line_start[0])*(line_end[1]-line_start[1])
    dotproduct = (point[0]-line_start[0])*(line_end[0]-line_start[0]) + (point[1]-line_start[1])*(line_end[1]-line_start[1])
    squaredlength = (line_end[0]-line_start[0])**2 + (line_end[1]-line_start[1])**2

    return (np.absolute(crossproduct)<=1e5)&(dotproduct>=0)&(dotproduct<=squaredlength)

# Computation
def TTC(samples, toreturn='dataframe'):
    if toreturn!='dataframe' and toreturn!='values':
        warnings.warn('Incorrect target to return. Please specify \'dataframe\' or \'values\'.')
    else:
        line_start = samples[['x_i','y_i']].values.T
        line_end = (samples[['x_i','y_i']].values + samples[['vx_i','vy_i']].values - samples[['vx_j','vy_j']].values).T

        direct_v = line_end - line_start
        perp_direct_v = np.array([-direct_v[1],direct_v[0]])
        direct_v_scale = np.tile(np.sqrt(direct_v[0]**2+direct_v[1]**2), (2,1))
        point_line00 = line_start + perp_direct_v/direct_v_scale*np.tile(samples.width_i.values, (2,1))/2
        point_line01 = point_line00 + direct_v
        point_line10 = line_start - perp_direct_v/direct_v_scale*np.tile(samples.width_i.values, (2,1))/2
        point_line11 = point_line10 + direct_v

        heading = samples[['hx_j','hy_j']].values
        perp_heading = np.array([-heading[:,1], heading[:,0]]).T
        heading_scale = np.tile(np.sqrt(heading[:,0]**2+heading[:,1]**2), (2,1)).T
        length = np.tile(samples.length_j.values, (2,1)).T
        width = np.tile(samples.width_j.values, (2,1)).T

        point_up = samples[['x_j','y_j']].values + heading/heading_scale*length/2
        point_down = samples[['x_j','y_j']].values - heading/heading_scale*length/2
        point1 = (point_up + perp_heading/heading_scale*width/2).T
        point2 = (point_up - perp_heading/heading_scale*width/2).T
        point3 = (point_down + perp_heading/heading_scale*width/2).T
        point4 = (point_down - perp_heading/heading_scale*width/2).T

        dist_mat = []
        ist_mat = []
        for point_a, point_b in zip([point1, point3, point1, point2],[point2, point4, point3, point4]):
            for point_line_a, point_line_b in zip([point_line00, point_line10],[point_line01, point_line11]):
                ### intersection point        
                ist = intersect(line(point_a, point_b), line(point_line_a, point_line_b))
                ist[:,~ison(point_a, point_b, ist)] = np.nan
                ist_mat.append(ist)
                ### distance from point to the relative velocity line
                dist = np.absolute((line_end[0]-line_start[0])*(line_start[1]-ist[1])-(line_start[0]-ist[0])*(line_end[1]-line_start[1]))/np.sqrt((line_end[0]-line_start[0])**2+(line_end[1]-line_start[1])**2)
                ### distance to start
                dist_ist = np.sqrt((ist[0]-line_start[0])**2+(ist[1]-line_start[1])**2-dist**2)
                dist_ist[np.isnan(dist_ist)] = np.inf
                dist_mat.append(dist_ist)

        dist2overlap = np.array(dist_mat).min(axis=0) - samples.length_i.values/2
        point2overlap = np.array(ist_mat)[np.array(dist_mat).argmin(axis=0), :, np.arange(len(samples))]
        direct_v = line_end - line_start
        direct_p = point2overlap.T - line_start
        leaving = (direct_v[0]*direct_p[0] + direct_v[1]*direct_p[1]) < 0
        dist2overlap[leaving] = np.inf

        TTC = dist2overlap/np.sqrt((samples.vx_i-samples.vx_j)**2+(samples.vy_i-samples.vy_j)**2)
        samples['TTC'] = TTC

        if toreturn=='dataframe':
            return samples
        elif toreturn=='values':
            return TTC