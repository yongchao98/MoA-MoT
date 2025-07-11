import numpy as np

def solve_packing_problem():
    """
    Solves the energy ball packing problem by determining the best container
    and the optimal number of balls that can be packed inside.
    """

    # Step 1: Find the optimal container.
    # Based on maximizing volume for a given surface area, the sphere is the best theoretical shape.
    # SA = 4 * pi * R^2 <= 1050  => R <= sqrt(1050 / (4*pi)) ~= 9.14 cm.
    # The largest radius that is a multiple of 0.5 cm is 9.0 cm.
    # Let's verify this choice.
    # Sphere(R=9.0): SA = 4*pi*9^2 = 1017.9 cm^2. Volume = 4/3*pi*9^3 = 3053.6 cm^3.
    # Cube(L=13.0): SA = 6*13^2 = 1014 cm^2. Volume = 13^3 = 2197 cm^3.
    # Cylinder(R=7.5, H=14.5): SA = 1036.5 cm^2. Volume = 2561.4 cm^3.
    # The sphere offers the most volume. We will proceed with the sphere of R=9.0 cm.
    
    container_description = "sphere r=9.0"
    container_radius = 9.0
    
    large_ball_radius = 2.0
    small_ball_radius = 1.0
    
    # Step 2 & 3: Pack large balls using a Simple Cubic (SC) lattice.
    # The lattice is centered at the origin (0,0,0). The spacing is determined
    # by the diameter of the large balls, d = 4.0 cm.
    # Center coordinates will be (4i, 4j, 4k) for integers i, j, k.
    # A large ball fits if its center `c` satisfies ||c|| + r_large <= R_container.
    # ||c|| <= 9.0 - 2.0 = 7.0.
    
    num_large_balls = 0
    large_ball_centers = []
    # We test integer indices i, j, k around the origin.
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                center = np.array([4.0 * i, 4.0 * j, 4.0 * k])
                if np.linalg.norm(center) <= 7.0:
                    num_large_balls += 1
                    large_ball_centers.append(center)
                    
    # The calculation shows that i,j,k in {-1, 0, 1} works, giving a 3x3x3 cube of balls.
    # For example, center (4,4,4) has norm sqrt(48) = 6.928 <= 7.0.
    # Center (8,0,0) has norm 8 > 7.0.
    # So, num_large_balls is 3*3*3 = 27.
    b = 27

    # Step 4: Pack small balls in the voids and near the boundary.
    num_small_balls = 0
    
    # a) Pack in the interstitial voids of the SC lattice.
    # These voids are centered at (2+4i, 2+4j, 2+4k).
    # A small ball fits if ||c|| + r_small <= R_container => ||c|| <= 8.0.
    # And it must not overlap with large balls. dist(c_small, c_large) >= r_small + r_large = 3.0.
    # The voids are at the center of 8 large balls.
    # e.g., void at (2,2,2) is sqrt(2^2+2^2+2^2)=3.46 from centers like (4,4,4), (0,0,0) etc.
    # 3.46 > 3.0, so it fits.
    # The SC lattice of 27 balls has 2x2x2 = 8 central voids.
    # Let's check if they are inside the sphere.
    # Void center (2,2,2) has norm 3.46. 3.46 + 1.0 <= 9.0. It fits.
    # All 8 voids at (+-2, +-2, +-2) fit.
    num_small_balls += 8
    
    # b) Pack in the space between the lattice and the spherical container wall.
    # Let's test positions on the axes, e.g., (x,0,0).
    # The outermost large ball center on the x-axis is at (4,0,0). It extends to x=6.
    # A small ball at (x_s,0,0) must be dist >= 3 from (4,0,0).
    # |x_s - 4| >= 3 => x_s >= 7 or x_s <= 1.
    # Let's test center (7,0,0). Norm is 7. 7 + r_small = 8 <= 9. It fits in the container.
    # Overlap with other large balls? Closest are (4,4,0) and (4,-4,0).
    # dist = sqrt((7-4)^2 + (4-0)^2) = sqrt(9+16) = 5 >= 3. It's valid.
    # We can place 6 such balls at (+-7,0,0), (0,+-7,0), (0,0,+-7).
    num_small_balls += 6
    
    a = num_small_balls
    
    # Step 5: Final result.
    print(f"Optimal container: {container_description}")
    print(f"Number of 1-cm balls (a): {a}")
    print(f"Number of 2-cm balls (b): {b}")
    print("\nFinal Answer Format:")
    
    # The prompt asks to output the numbers in the final equation.
    # Although this is unusual, it likely refers to the energy calculation.
    large_ball_energy = 10
    small_ball_energy = 1
    total_energy = a * small_ball_energy + b * large_ball_energy
    
    final_answer_str = f"[{container_description}]{a};{b}"
    
    print(f"The final configuration is '{final_answer_str}'")
    print(f"The total energy is calculated as: {a} * {small_ball_energy} MJ + {b} * {large_ball_energy} MJ = {total_energy} MJ.")

    print(f"\n<<<[{container_description}]{a};{b}>>>")

solve_packing_problem()
