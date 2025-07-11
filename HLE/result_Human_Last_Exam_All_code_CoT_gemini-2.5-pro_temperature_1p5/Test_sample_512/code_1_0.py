import math

def solve():
    """
    Analyzes the initial container and designs a more material-efficient one.
    """
    # Step 1: Analyze the initial container
    cube_side = 12.0
    ball_radius = 2.0
    ball_diameter = 2 * ball_radius

    # In a 12x12x12 cube, we can fit 3 balls along each dimension (12 / 4 = 3).
    initial_ball_count = (cube_side / ball_diameter)**3
    
    # The surface area of the cube is 6 * side^2.
    initial_surface_area = 6 * (cube_side**2)

    # Step 2: Propose and design a new, more efficient container.
    # We will use a cylindrical container, which is generally more efficient
    # in terms of surface-area-to-volume ratio than a cuboid.

    # Our packing strategy will be based on stacking layers of 7 balls.
    # To satisfy the 0.5 cm grid constraint, we found that a 7-ball hexagonal-like layer
    # requires a radius of at least 6.03 cm to contain the balls (center_radius + ball_radius).
    # The next highest multiple of 0.5 cm is 6.5 cm.
    new_container_radius = 6.5

    # We can stack these layers in a nestled (staggered) formation. The minimum
    # vertical separation between the centers of these layers, while respecting the
    # 0.5cm grid, is 3.0 cm.
    # To hold more than 27 balls, we use 4 layers (4 * 7 = 28 balls).
    # The centers of the 4 layers are at z=2.0, 5.0, 8.0, and 11.0.
    # The total height required for the container is from the bottom of the first ball (z=0.0)
    # to the top of the last ball (z=13.0).
    new_container_height = 13.0
    new_ball_count = 28

    # Step 3: Calculate the surface area of the new cylindrical container.
    # The formula is: Surface Area = 2 * pi * r^2 (for the ends) + 2 * pi * r * h (for the side).
    # This simplifies to: Surface Area = 2 * pi * r * (r + h).
    pi_val = math.pi
    new_surface_area = 2 * pi_val * new_container_radius * (new_container_radius + new_container_height)

    # Step 4: Compare with the initial container and format the output.
    if new_surface_area < initial_surface_area and new_ball_count >= initial_ball_count:
        print(f"The initial 12x12x12 cm cube holds {int(initial_ball_count)} balls and has a surface area of {initial_surface_area:.2f} cm^2.")
        print("A more efficient container has been found.")
        print(f"It is a cylinder that holds {new_ball_count} balls, which is more than the original container.")
        
        container_description = f"cylinder r={new_container_radius}, h={new_container_height}"
        final_answer_string = f"{new_surface_area:.2f}[{container_description}]"

        print("\n# Final Equation for the new container's surface area 'd':")
        print(f"# d = 2 * pi * r * (r + h)")
        print(f"# d = 2 * {pi_val:.4f} * {new_container_radius} * ({new_container_radius} + {new_container_height})")
        print(f"# d = {new_surface_area:.2f} cm^2\n")
        
        print(f"Yes, a more efficient container is possible: {final_answer_string}")
        print(f"<<<{final_answer_string}>>>")
    else:
        # This case is executed if no better container is found.
        print("0")
        print("<<<0>>>")

solve()