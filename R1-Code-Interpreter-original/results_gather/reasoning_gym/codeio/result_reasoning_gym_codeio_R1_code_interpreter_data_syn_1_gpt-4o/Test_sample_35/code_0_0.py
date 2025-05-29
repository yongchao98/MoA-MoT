import math
import numpy as np

class Obstacle:
    def __init__(self, x=0, y=0, r=0):
        self.x = x
        self.y = y
        self.r = abs(r)

    def does_line_segment_intersect(self, x1=0, y1=0, x2=0, y2=0):
        point_a = np.array([x1, y1])
        point_b = np.array([x2, y2])
        point_circle = np.array([self.x, self.y])
        direction_vector = point_b - point_a
        circle_out_vector = point_a - point_circle

        A = np.dot(direction_vector, direction_vector)
        B = np.dot(circle_out_vector, direction_vector) * 2
        C = np.dot(circle_out_vector, circle_out_vector) - pow(self.r, 2)

        discrim = pow(B, 2) - (4 * A * C)

        if discrim < 0:
            return False

        discrim = math.sqrt(discrim)
        t1 = (-B - discrim) / (2 * A)
        t2 = (-B + discrim) / (2 * A)
        if (t1 >= 0 and t1 <= 1) or (t2 >= 0 and t2 <= 1):
            return True

        return False

    def is_at_least_one_endpoint_within_radius(self, x1=0, y1=0, x2=0, y2=0):
        d1 = math.sqrt((self.x - x1) * (self.x - x1) + (self.y - y1) * (self.y - y1))
        d2 = math.sqrt((self.x - x2) * (self.x - x2) + (self.y - y2) * (self.y - y2))
        return d1 <= self.r or d2 <= self.r

    def does_obstacle_envelope_edge(self, x1=0, y1=0, x2=0, y2=0):
        d1 = math.sqrt((self.x - x1) * (self.x - x1) + (self.y - y1) * (self.y - y1))
        d2 = math.sqrt((self.x - x2) * (self.x - x2) + (self.y - y2) * (self.y - y2))
        return d1 <= self.r and d2 <= self.r

    def is_obstacle_on_edge(self, x1=0, y1=0, x2=0, y2=0):
        does_segment_intersect = self.does_line_segment_intersect(x1, y1, x2, y2)
        is_line_enveloped = self.does_obstacle_envelope_edge(x1, y1, x2, y2)
        at_least_once = self.is_at_least_one_endpoint_within_radius(x1, y1, x2, y2)
        return does_segment_intersect or is_line_enveloped or at_least_once

def main_solution(circle_x, circle_y, circle_r, line_x1, line_y1, line_x2, line_y2):
    obstacle = Obstacle(circle_x, circle_y, circle_r)
    result = obstacle.is_obstacle_on_edge(line_x1, line_y1, line_x2, line_y2)
    return {"intersects": result}

# Test the function with the proposed input
result = main_solution(0, 0, 1, 2, 2, 3, 3)
print(result)