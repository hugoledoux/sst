import os
import sys
import time
import psutil

import startinpy

import numpy as np

from multiprocessing import cpu_count, Process, Lock, Queue, current_process
from scipy.spatial import KDTree

COARSE_THRESHOLD = 2
FINE_THRESHOLD = 0.2


class MemoryUsage:
    def __init__(self, process_name, timestamp, memory_usage):
        self.process_name = process_name
        self.timestamp = timestamp
        self.memory_usage = memory_usage


class Triangulation:
    def __init__(self):
        self.cell_size = None

        self.min_x = None
        self.min_y = None
        self.max_x = None
        self.max_y = None

    def set_bbox(self, min_x, min_y, max_x, max_y):
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = max_x
        self.max_y = max_y

    def finalize(self, input_line, grid_x, grid_y, vertices, lock, memory_usage_queue):
        stdout_lines = []
        if len(vertices) > 0:

            triangulation = startinpy.DT()

            x_vals = []
            y_vals = []
            z_vals = []

            for vertex_id, point in vertices.items():
                x_vals.append(point[0])
                y_vals.append(point[1])
                z_vals.append(point[2])

            tree = KDTree(np.c_[x_vals, y_vals])

            corner_points = [
                [self.min_x + (self.cell_size * grid_x), self.min_y + (self.cell_size * grid_y)],
                [self.min_x + (self.cell_size * grid_x) + self.cell_size - 1E-5, self.min_y + (self.cell_size * grid_y)],
                [self.min_x + (self.cell_size * grid_x), self.min_y + (self.cell_size * grid_y) + self.cell_size - 1E-5],
                [self.min_x + (self.cell_size * grid_x) + self.cell_size - 1E-5, self.min_y + (self.cell_size * grid_y) + self.cell_size - 1E-5]
            ]

            near_corner_points = []

            for corner_point in corner_points:
                # Get nearest point to corner
                distances, indexes = tree.query(corner_point, k=10)

                queried_z_vals = [z_vals[index] for index in indexes if index < len(z_vals)]

                # add a corner point with average z value of 10 nearest
                near_corner_points.append([corner_point[0], corner_point[1], sum(queried_z_vals) / len(queried_z_vals)])

            triangulation.insert(near_corner_points)

            to_delete = []

            # First coarse loop

            for key, vertex in vertices.items():
                x = vertex[0]
                y = vertex[1]
                z = vertex[2]

                try:
                    interpolated_value = triangulation.interpolate_tin_linear(x, y)

                    if abs(interpolated_value - z) > COARSE_THRESHOLD:
                        triangulation.insert_one_pt(x, y, z)
                        to_delete.append(key)

                # In rare cases we get a point outside CH due to ----00.0 being counted as wrong cell
                # FIXME: Adjust get_cell function to return correct cell for ----00.0 points
                except OSError:
                    pass

            for key in reversed(to_delete):
                del vertices[key]

            # Fine loop

            for key, vertex in vertices.items():
                x = vertex[0]
                y = vertex[1]
                z = vertex[2]

                try:
                    interpolated_value = triangulation.interpolate_tin_linear(x, y)

                    if abs(interpolated_value - z) > FINE_THRESHOLD:
                        triangulation.insert_one_pt(x, y, z)

                # In rare cases we get a point outside CH due to ----00.0 being counted as wrong cell
                # FIXME: Adjust get_cell function to return correct cell for ----00.0 points
                except OSError:
                    pass

            if triangulation.number_of_vertices() > 4:
                for i in [1, 2, 3, 4]:
                    triangulation.remove(i)

            for vertex in triangulation.points:
                # Don't print infinite vertex
                if vertex[0] > 0:
                    stdout_lines.append("v " + str(vertex[0]) + " " + str(vertex[1]) + " " + str(vertex[2]) + "\n")

        memory_usage_queue.put(MemoryUsage(current_process().name, round(time.time()), psutil.Process(os.getpid()).memory_info().rss))

        stdout_lines.append(input_line)

        with lock:
            sys.stdout.write("".join(stdout_lines))
            sys.stdout.flush()

        sys.stderr.write(current_process().name + " - FINISHED.\n")

class Processor:
    def __init__(self, dt):
        self.triangulation = dt

        self.vertex_id = 1
        self.vertices = {}

        # self.sprinkling = False
        self.processes = []

        self.last_log_time = round(time.time())

        self.stdout_lock = Lock()
        self.memory_usage_queue = Queue()

        self.memory_usage_queue.put(MemoryUsage("Main", self.last_log_time, psutil.Process(os.getpid()).memory_info().rss))

        self.memory_usage_writer = Process(target=self.write_memory_usage, args=(self.memory_usage_queue,), daemon=True)
        self.memory_usage_writer.start()

    def write_memory_usage(self, memory_usage_queue):
        while True:
            with open(os.path.join(os.getcwd(), "memlog_direct_refinement.csv"), "a") as memory_log_file:
                val = memory_usage_queue.get()

                if val:
                    memory_log_file.write(str(val.process_name) + ", " + str(val.timestamp) + ", " + str(val.memory_usage) + "\n")
                else:
                    time.sleep(0.5)

    def process_line(self, input_line):
        split_line = input_line.rstrip("\n").split(" ")

        identifier = split_line[0]
        data = split_line[1:]

        current_time = round(time.time())

        if current_time != self.last_log_time:
            self.memory_usage_queue.put(MemoryUsage("Main", current_time, psutil.Process(os.getpid()).memory_info().rss))
            self.last_log_time = current_time

        if identifier == "#" or identifier == "":
            pass

        elif identifier == "n":
            # Total number of points
            self.triangulation.total_points = int(data[0])
            sys.stdout.write(input_line)

        elif identifier == "c":
            # Grid dimensions (cXc)
            sys.stdout.write(input_line)

        elif identifier == "s":
            # Cell size
            self.triangulation.cell_size = int(data[0])
            sys.stdout.write(input_line)

        elif identifier == "b":
            # bbox
            self.triangulation.set_bbox(float(data[0]), float(data[1]), float(data[2]), float(data[3]))
            sys.stdout.write(input_line)

            sys.stderr.write(input_line)
            sys.stderr.flush()

        elif identifier == "v":
            self.vertices[self.vertex_id] = [float(data[0]), float(data[1]), float(data[2])]
            self.vertex_id += 1

        elif identifier == "w":
            sys.stdout.write(input_line)

        elif identifier == "x":
            # cell finalizer
            # While sprinkling, don't bother processing since all finalized cells now are still empty anyways
            # if self.sprinkling:
            #     sys.stdout.write(input_line)
            #     return

            sys.stderr.write("Starting new process to finalize cell: {}, {}. Processing currently running: {}\n".format(data[0], data[1], len(self.processes)))
            sys.stderr.flush()

            sleep_time = 1

            # Ensure total number of processes never exceeds capacity
            while len(self.processes) >= cpu_count() - 2:
                for i in reversed(range(len(self.processes))):
                    if not self.processes[i].is_alive():
                        del self.processes[i]

                time.sleep(sleep_time)

            process = Process(target=self.triangulation.finalize, args=(input_line, int(data[0]), int(data[1]), self.vertices, self.stdout_lock, self.memory_usage_queue,), daemon=True)
            self.vertices = {}
            self.vertex_id = 1
            self.processes.append(process)
            process.start()

        else:
            # Unknown identifier in stream
            pass

        sys.stdout.flush()


if __name__ == "__main__":
    triangulation = Triangulation()
    processor = Processor(triangulation)

    start_time = time.time()

    for stdin_line in sys.stdin:
        processor.process_line(stdin_line)

    for process in processor.processes:
        process.join()

    # sys.stderr.write("duration: " + str(time.time() - start_time) + "\n")
