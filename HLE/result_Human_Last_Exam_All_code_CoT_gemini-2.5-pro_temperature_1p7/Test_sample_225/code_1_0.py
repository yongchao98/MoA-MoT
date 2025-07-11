import numpy as np
import math

def find_intersections():
    """
    Finds the number of self-intersection points of the Lissajous curve
    x = cos(a*t), y = sin(b*t) using a numerical simulation.
    """
    # The parameters for the curve (cos(a*t), sin(b*t))
    a = 9
    b = 5

    # Simulation parameters
    grid_resolution = 4000  # How fine the grid is
    time_steps = 200000  # Number of steps to trace the curve
    t_min_diff = 0.1  # Minimum time difference to be considered a new crossing

    # Data structures
    # A dictionary to store the time of the first visit to a grid cell
    grid = {}
    # A set to store the coordinates of unique intersection points
    intersection_points = set()

    # The curve's full period is 2*pi since gcd(a,b)=1
    t_values = np.linspace(0, 2 * np.pi, time_steps, endpoint=False)

    for t in t_values:
        x = math.cos(a * t)
        y = math.sin(b * t)

        # Convert continuous coordinates [-1, 1] to discrete grid coordinates
        # N/2 is a scaling factor, adding N/2 shifts the origin
        ix = int((grid_resolution / 2) * (x + 1))
        iy = int((grid_resolution / 2) * (y + 1))
        
        # Check if we've visited this grid cell before
        if (ix, iy) in grid:
            t_old = grid[(ix, iy)]
            # If the time difference is large enough, it's a real self-intersection
            if abs(t - t_old) > t_min_diff:
                intersection_points.add((ix, iy))
        else:
            # First time visiting this cell, record the time
            grid[(ix, iy)] = t

    num_intersections = len(intersection_points)
    
    # The final equation is not a formula, but a statement of the result
    # We output the numbers 'a' and 'b' and the final count as requested.
    print(f"For the curve defined by x = cos({a}t) and y = sin({b}t):")
    print(f"The number of self-intersection points is {num_intersections}")


find_intersections()
<<<32>>>