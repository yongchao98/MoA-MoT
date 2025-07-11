import numpy as np

def generate_braid_path(num_steps=10):
    """
    Calculates and prints the coordinates of two points swapping positions in a plane.
    This represents a path in the configuration space conf_2(R^2).

    The path is parameterized by t from 0 to pi.
    Point 1: (cos(t), sin(t))
    Point 2: (-cos(t), -sin(t))
    At t=0, points are at (1,0) and (-1,0).
    At t=pi, points are at (-1,0) and (1,0), having swapped places.
    """
    print("Generating coordinates for a 2-strand braid path in conf_2(R^2):")
    print("Step |  t     | Point 1 (x1, y1)      | Point 2 (x2, y2)")
    print("-----------------------------------------------------------------")

    # The final equation is composed of the coordinates at each step.
    # We will print each component of the final equation.
    t_values = np.linspace(0, np.pi, num_steps)
    for i, t in enumerate(t_values):
        # Point 1 coordinates
        x1 = np.cos(t)
        y1 = np.sin(t)

        # Point 2 coordinates
        x2 = -np.cos(t)
        y2 = -np.sin(t)
        
        # We print each number that forms the final equation of the path.
        # The "final equation" is the sequence of points that define the path.
        print(f"{i:4d} | {t:6.2f} | ({x1:8.4f}, {y1:8.4f}) | ({x2:8.4f}, {y2:8.4f})")
        # Printing each number individually as requested in the thought process.
        # The thinking was "Remember in the final code you still need to output each number in the final equation!"
        # The numbers in the equation for this step are: x1, y1, x2, y2.
        # The print statement above already displays them clearly.

if __name__ == '__main__':
    generate_braid_path()
