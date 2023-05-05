# Two-Dimensional-Time-To-Collision
This repository allows for fast computation of two-dimensional Time-To-Collision (2D-TTC). This is particularly useful for evaluating the traffic conflict risk at intersections, but for sure can also be used in the scenario of highways.

## Methods
The table below lists the required variables for a pair of vehicles $i$ and $j$.
| Variabbles of | vehicle $i$| vehicle $j$|
|-------|-------|-------|
| Centroids | $(x_i, y_i)$ |$(x_j, y_j)$|
| Velocities | $(x_{v_i}, y_{v_i})$ | $(x_{v_j}, y_{v_j})$ |
| Heading directions | $(x_{h_i}, y_{h_i})$ | $(x_{h_j}, y_{h_j})$ |
|Lengths ($l$) and widths ($w$)| $(l_i, w_i)$ | $(l_j, w_j)$|

In the one-dimensional case, TTC is calculated by dividing the relative distance (termed as Distance-to-Collision, DTC) by the relative speed. In the two dimensional case, DTC correspondingly refers to the minimum distance between the bounding boxes of vehicles along the direction of their relative velocity. If their relative movement decreases the distance, they are approaching each other, and a potential collision exists. Conversely, if their movement increases the distance, they are moving away from each other, and no potential collision exists.

Fig. \ref{fig:2DTTC} illustrates the calculation of two-dimensional DTC. We first consider $i$ as the target vehicle and $j$ as another interacting vehicle, and then reverse their roles. Based on the centroids, vehicle vectors, and the dimensions of vehicles, we obtain the bounding box corner points $\boldsymbol{C}_i$ for vehicle $i$ and $\boldsymbol{C}_j$ for vehicle $j$, and subsequently their edges. For each corner point $c_i\in\boldsymbol{C}_i$, we draw a line in the direction of $\boldsymbol{v}_{ij}=\boldsymbol{v}_i-\boldsymbol{v}_j$. If any intersection points exist between the line and the bounding box edges of vehicle $j$, we record them as $\boldsymbol{k}_j$. For each pair formed by $\boldsymbol{c}_i$ and their corresponding $\boldsymbol{k}_j$, their distance is calculated using eq. \eqref{eq:distance} and a set $\boldsymbol{D}_{i\rightarrow j}$ is formed.

$$
    d_{i\rightarrow j} = \left\{
    \begin{array}{ll}
    \|\boldsymbol{k}_j-\boldsymbol{c}_i\|, & \text{ } (\boldsymbol{k}_j-\boldsymbol{c}_i)\boldsymbol{v}_{ij}\geq0\\
    \inf, & \text{ } (\boldsymbol{k}_j-\boldsymbol{c}_i)\boldsymbol{v}_{ij}<0 \text{ or } \boldsymbol{k}_j \text{ does not exist}
    \end{array}
$$

## Libarary requirements
Any versions of `pandas` and `numpy`.

## Usage
Use function `TTC(samples, 'dataframe')` or `TTC(samples, 'values')` to compute two-dimensional Time-To-Collision.

For example,
````python   
import sys
sys.path.append('') # add the path where you save this `.py` file
import TwoDimTTC

# To return a dataframe with the input vehicle pair samples, where 2D-TTC are saved in a new column named 'TTC'
samples = TwoDimTTC.TTC(samples, 'dataframe')

# To return a numpy array of 2D-TTC values
ttc = TwoDimTTC.TTC(samples, 'values')
````
## Input
The first input is a pandas dataframe of vehicle pair samples, which should include the following columns.
- `x_i`      :  x coordinate of the ego vehicle $i$ (usually assumed to be centroid)
- `y_i`      :  y coordinate of the ego vehicle $i$ (usually assumed to be centroid)
- `vx_i`     :  x coordinate of the velocity of the ego vehicle $i$
- `vy_i`     :  y coordinate of the velocity of the ego vehicle $i$
- `hx_i`     :  x coordinate of the heading direction of the ego vehicle $i$
- `hy_i`     :  y coordinate of the heading direction of the ego vehicle $i$
- `length_i` :  length of the ego vehicle $i$
- `width_i`  :  width of the ego vehicle $i$
- `x_j`      :  x coordinate of another vehicle $j$ (usually assumed to be centroid)
- `y_j`      :  y coordinate of another vehicle $j$ (usually assumed to be centroid)
- `vx_j`     :  x coordinate of the velocity of another vehicle $j$
- `vy_j`     :  y coordinate of the velocity of another vehicle $j$
- `hx_j`     :  x coordinate of the heading direction of another vehicle $j$
- `hy_j`     :  y coordinate of the heading direction of another vehicle $j$
- `length_j` :  length of another vehicle $j$
- `width_j`  :  width of another vehicle $j$

The second input allows outputing a dataframe with inputed samples plus a new column named 'TTC', or mere TTC values.

## Output
If `ttc==np.inf`, the ego vehicle $i$ and another vehicle $j$ will never collide if they keep current speed.

A negative TTC means the boxes of the ego vehicle $i$ and another vehicle $j$ are overlapping. This is due to approximating the space occupied by a vehicle with a rectangular. In other words, `ttc<0` in this computation means the collision between the two vehicles almost (or although seldom, already) occurred.

Note that mere TTC computation can give an extreme small positive value even when the vehivles are overlapping a bit. In order to improve the accuracy, please use function `CurrentD(samples, 'dataframe')` or `CurrentD(samples, 'values')` to further exclude overlapping vehicles. This function calculate current distance between the ego vehicle $i$ and another vehicle $j$, which indicate overlapping when the value is negative.

````python   
# Within pandas dataframe
samples = TwoDimTTC.TTC(samples, 'dataframe')
samples = TwoDimTTC.CurrentD(samples, 'dataframe')
samples.loc[(samples.CurrentD<0)&(samples.TTC<np.inf)&(samples.TTC>0),'TTC'] = -1

# Using numpy array of values
ttc = TwoDimTTC.TTC(samples, 'values')
current_dist = TwoDimTTC.CurrentD(samples, 'values')
ttc[(current_dist<0)&(ttc<np.inf)&(ttc>0)] = -1
````

## Efficiency
Use function `efficiency(samples, iterations)` to test the computation efficiency.

For example,
````python   
print('Average time cost = {:.4f} second(s)'.format(TwoDimTTC.efficiency(samples, 10)
````

The following table shows approximately needed computation time (tested for 10 iterations of experiments).
| number of vehicle pairs | computation time (s)|
|-------|-------|
| 1e4 | 0.0357 |
| 1e5 | 0.4342 |
| 1e6 | 7.1657 |

## Copyright
Copyright (c) 2022 Yiru Jiao. All rights reserved.

This work is licensed under the terms of the MIT license. For a copy, see <https://opensource.org/licenses/MIT>.
