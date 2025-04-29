# Given a file with the following contents:
# Each line contains multiple points separated by space.
# The points are in the format (x,y).
# Each line is one partition.
# Plot all the points in the same figure and draw boundaries for each partition.
#
# Example input:
# (0.0898955,0.33079) (0.685802,0.903268) 
# (0.0540748,0.923177) (0.0740991,0.928702) 
# (0.203431,0.358879) (0.206965,0.699515) 
# (0.722753,0.629212) (0.720795,0.730465) 
# (0.884282,0.0917934) (0.884531,0.463237)

import matplotlib.pyplot as plt
import re
import sys

def main():
  colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
  pt_size = 2
  with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
    colors_idx = 0
    for line in lines:
      # if the line does not start with '(' or is empty, skip it
      if not line.strip() or not line.startswith('('):
        continue
      points = line.rstrip().split(' ')
      points_list = []
      maxx = maxy = -1
      minx = miny = 1000000
      for pt in points:
        pt = pt.strip('()')
        x, y = pt.split(',')
        points_list.append((float(x), float(y)))
        maxx = max(maxx, float(x))
        maxy = max(maxy, float(y))
        minx = min(minx, float(x))
        miny = min(miny, float(y))
      c = colors[colors_idx % len(colors)]
      colors_idx += 1
      # Draw the points and make the points size smaller
      plt.scatter(*zip(*points_list), color=c, s=pt_size)
      # Draw the boundaries
      plt.plot([minx, maxx], [miny, miny], c)
      plt.plot([minx, maxx], [maxy, maxy], c)
      plt.plot([minx, minx], [miny, maxy], c)
      plt.plot([maxx, maxx], [miny, maxy], c)
      
      plt.savefig(sys.argv[2])

if __name__ == '__main__':
  main()