import csv

from matplotlib import pyplot as plt


def get_stats(filename):
  with open(filename, "r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    num_columns = 9
    stats = [[] for i in range(9)]
    for row in csv_reader:
      for i in range(num_columns):
        stats[i].append(float(row[i]))
    return stats


def main():
  stats_lipp = get_stats('lipp.csv')
  stats_mlipp = get_stats('mlipp.csv')
  stats_kdtree = get_stats('kdtree.csv')
  stats_zorder = get_stats('zorder.csv')

  graph_titles = [
    "Construction (Insert)", 
    "Point Query (Insert)",
    "Max Depth (Insert)",
    "Avg Depth (Insert)",
    "Construction (Bulk)", 
    "Point Query (Bulk)",
    "Max Depth (Bulk)",
    "Avg Depth (Bulk)"
  ]

  for i in range(len(graph_titles)):
    plt.plot(stats_lipp[0], stats_lipp[i+1], label="lipp")
    plt.plot(stats_mlipp[0], stats_mlipp[i+1], label="mlipp")
    plt.plot(stats_kdtree[0], stats_kdtree[i+1], label="kdtree")
    plt.plot(stats_zorder[0], stats_zorder[i+1], label="zorder")
    plt.ylim(0, None)
    plt.legend()
    plt.title(graph_titles[i])
    plt.show()

if __name__ == '__main__':
  main()