import numpy as np
import matplotlib.pyplot as plt

file_path = r"./computation_result.txt"
# file_path = r"./build/computation_result.txt"

N_x = 400
N_y = 600
x_min = -0.1
x_max = 1.1
y_min = -1.1
y_max = 1.1

def reader(file_path):
    data = np.empty((N_x + 1, N_y + 1))
    with open(file_path) as file:
        time, iteration_no = file.readline().split()
        time = float(time)
        iteration_no = int(iteration_no)
        for i in range(N_x + 1):
            data[i, :] = list(map(float, file.readline().split()))
    return time, iteration_no, data

def plot_solution(data):
    _, ax = plt.subplots()
    ax.imshow(data.T, cmap = "RdYlBu", extent = (x_min, x_max, y_min, y_max), origin = "lower")
    ax.set_aspect('equal')
    plt.show()

if __name__ == "__main__":
    time, iteration_no, data = reader(file_path)
    plot_solution(data)
