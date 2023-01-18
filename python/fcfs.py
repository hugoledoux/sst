import os
import sys
import time
import psutil

import startinpy

import numpy as np

from scipy.spatial import KDTree

COARSE_THRESHOLD = 2
FINE_THRESHOLD = 0.2


# class MemoryUsage:
#     def __init__(self, process_name, timestamp, memory_usage):
#         self.process_name = process_name
#         self.timestamp = timestamp
#         self.memory_usage = memory_usage


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

    def finalize(self, input_line, grid_x, grid_y, vertices):
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

        # memory_usage_queue.put(MemoryUsage(current_process().name, round(time.time()), psutil.Process(os.getpid()).memory_info().rss))

        stdout_lines.append(input_line)

        sys.stdout.write("".join(stdout_lines))
        sys.stdout.flush()

        # sys.stderr.write(current_process().name + " - FINISHED.\n")




if __name__ == "__main__":
    triangulation = Triangulation()
    vertices = {}
    vertex_id = 1

    for input_line in sys.stdin:
        split_line = input_line.rstrip("\n").split(" ")

        identifier = split_line[0]
        data = split_line[1:]

        if identifier == "#" or identifier == "":
            pass

        elif identifier == "n":
            # Total number of points
            triangulation.total_points = int(data[0])
            sys.stdout.write(input_line)

        elif identifier == "c":
            # Grid dimensions (cXc)
            sys.stdout.write(input_line)

        elif identifier == "s":
            # Cell size
            triangulation.cell_size = int(data[0])
            sys.stdout.write(input_line)

        elif identifier == "b":
            # bbox
            triangulation.set_bbox(float(data[0]), float(data[1]), float(data[2]), float(data[3]))
            sys.stdout.write(input_line)

            sys.stderr.write(input_line)
            sys.stderr.flush()

        elif identifier == "v":
            vertices[vertex_id] = [float(data[0]), float(data[1]), float(data[2])]
            vertex_id += 1

        elif identifier == "w":
            sys.stdout.write(input_line)

        elif identifier == "x":

            triangulation.finalize(input_line, int(data[0]), int(data[1]), vertices)
            vertices = {}
            vertex_id = 1
            # self.processes.append(process)
            # process.start()

        else:
            # Unknown identifier in stream
            pass

        sys.stdout.flush()

