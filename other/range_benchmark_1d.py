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
    lipp = [[], []]
    kdtree = [[], []]
    for row in csv_reader:
      if row[0] == "lipp":
        lipp[0].append(int(row[1]))
        lipp[1].append(float(row[2]))
      else: # row[0] == "kdtree"
        kdtree[0].append(int(row[1]))
        kdtree[1].append(float(row[2]))
    return group_by(order_by(lipp)), group_by(order_by(kdtree))

def main():
  lipp, kdtree = get_stats('range_benchmark_1d.csv')

  print(lipp)

  plt.plot(lipp[0], lipp[1], label="lipp")
  plt.plot(kdtree[0], kdtree[1], label="kdtree")
  plt.ylim(0, None)
  plt.legend()
  plt.title("Range query duration")
  plt.show()

if __name__ == '__main__':
  main()