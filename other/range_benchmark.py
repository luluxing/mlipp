import csv

from matplotlib import pyplot as plt


def round_to(x, base):
    return base * round(x/base)

def group_by(index):
  # return index

  xs = index[0]
  ys = index[1]

  bucket = 1

  new_xs = [round_to(xs[0], bucket)]
  new_ys = [ys[0]]
  n = 1
  for i in range(1, len(xs)):
    if (round_to(xs[i], bucket) <= round_to(xs[i-1], bucket)):
      new_ys[-1] += ys[i]
      n += 1
    else:
      new_ys[-1] /= n
      new_xs.append(round_to(xs[i], bucket))
      new_ys.append(ys[i])
      n = 1
  new_ys[-1] /= n
  return [new_xs, new_ys]

def order_by(index):
  xs = index[0]
  ys = index[1]

  n = len(xs)
  xy = [(xs[i], ys[i]) for i in range(n)]
  xy = sorted(xy)

  xs = [xy[i][0] for i in range(n)]
  ys = [xy[i][1] for i in range(n)]
  return [xs, ys]

def get_stats(filename):
  with open(filename, "r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    mlipp_kd = [[], []]
    mlipp_zorder = [[], []]
    kdtree = [[], []]
    for row in csv_reader:
      if row[0] == "mlipp_kd":
        mlipp_kd[0].append(int(row[1]))
        mlipp_kd[1].append(float(row[2]))
      elif row[0] == "mlipp_zorder":
        mlipp_zorder[0].append(int(row[1]))
        mlipp_zorder[1].append(float(row[2]))
      else: # row[0] == "kdtree"
        kdtree[0].append(int(row[1]))
        kdtree[1].append(float(row[2]))
    return order_by(group_by(mlipp_kd)), order_by(group_by(mlipp_zorder)), order_by(group_by(kdtree))

def main():
  mlipp_kd, mlipp_zorder, kdtree = get_stats('range_benchmark_medium.csv')

  plt.plot(mlipp_kd[0], mlipp_kd[1], label="mlipp_kd")
  plt.plot(mlipp_zorder[0], mlipp_zorder[1], label="mlipp_zorder")
  plt.plot(kdtree[0], kdtree[1], label="kdtree")
  plt.ylim(0, None)
  plt.legend()
  plt.title("Range query duration")
  plt.show()

if __name__ == '__main__':
  main()