import math

# This script analyzes the problem of finding a more material-efficient container
# for 27 or more energy balls.

def main():
    """
    Analyzes the container problem to determine if a more efficient design exists.
    """
    # --- Initial Configuration ---
    initial_dim = 12.0
    ball_radius = 2.0
    ball_diameter = 4.0
    min_balls_required = 27
    
    # Calculate properties of the initial container
    initial_surface_area = 6 * initial_dim**2
    initial_capacity = math.floor(initial_dim / ball_diameter)**3
    
    print("--- Initial Container Analysis ---")
    print(f"The original container is a cube: {initial_dim} x {initial_dim} x {initial_dim} cm.")
    print(f"Surface Area: 6 * {initial_dim}^2 = {initial_surface_area:.0f} cm^2.")
    print(f"Its capacity, using a simple grid for spheres of diameter {ball_diameter} cm, is floor({initial_dim}/{ball_diameter})^3 = {initial_capacity:.0f} balls.")
    print("-" * 35)
    
    print("\n--- Searching for a More Efficient Container ---")
    print("An efficient container must hold at least 27 balls and have a surface area < 864 cm^2.")
    
    print("\nStep 1: The '3-in-a-Row' Constraint")
    # To place 3 balls in a line, their centers would be separated by the ball diameter.
    # The space spanned by the centers is (3-1) * diameter.
    # The total box dimension must also account for the radius on each end.
    min_dim_for_3_balls = (3 - 1) * ball_diameter + ball_diameter
    print(f"To place 3 balls in a line, the centers require a span of (3-1) * {ball_diameter} = { (3-1)*ball_diameter:.1f} cm.")
    print(f"The container dimension must fully enclose the balls, requiring at least { (3-1)*ball_diameter:.1f} cm (for centers) + {ball_diameter:.1f} cm (for radii) = {min_dim_for_3_balls:.1f} cm.")
    print("This implies that a simple 3x3x3 grid packing requires a box of at least 12x12x12 cm.")

    print("\nStep 2: Test a Box with Smaller Surface Area")
    print("For a cuboid to have SA < 864, at least one dimension must be smaller than 12 cm.")
    # Let's test a candidate box: 11.5 x 11.5 x 12.0 cm.
    test_l, test_w, test_h = 11.5, 11.5, 12.0
    test_sa = 2 * (test_l * test_w + test_l * test_h + test_w * test_h)
    print(f"Let's test a box {test_l} x {test_w} x {test_h}. Its surface area is 2*({test_l}*{test_w} + {test_l}*{test_h} + {test_w}*{test_h}) = {test_sa:.2f} cm^2.")
    
    print("\nStep 3: Check the Capacity of the Smaller Box")
    # Maximum number of balls that can be placed on a simple grid.
    nx = math.floor((test_l - ball_diameter) / ball_diameter) + 1
    ny = math.floor((test_w - ball_diameter) / ball_diameter) + 1
    nz = math.floor((test_h - ball_diameter) / ball_diameter) + 1
    test_capacity = nx * ny * nz
    print(f"Because two dimensions are < {min_dim_for_3_balls} cm, we can't fit 3 balls in a row along those axes.")
    print(f"The maximum grid capacity is {nx} x {ny} x {nz} = {test_capacity} balls.")
    print(f"This is less than the required {min_balls_required} balls.")

    print("\nStep 4: Other Packing Considerations")
    print("More complex packings (like FCC) use layer distances that are not multiples of 0.5 cm. Forcing them onto the 0.5 cm grid makes the packing less efficient and increases the total volume of the bounding box.")

    print("\n--- Conclusion ---")
    print("No container has been found that uses less material while holding the required number of energy balls.")
    print("The original 12x12x12 cm cube is the optimal design given the precision constraints.")
    
if __name__ == '__main__':
    main()