def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """
    
    # Step 1: Define the Betti numbers for the 2-torus (T^2).
    # The Betti numbers describe the topology of the space.
    # b_0: number of connected components.
    # b_1: number of 1-dimensional "holes" or "tunnels".
    # b_2: number of 2-dimensional "voids".
    b0 = 1  # The torus is one connected component.
    b1 = 2  # The torus has two independent circular loops.
    b2 = 1  # The torus encloses one 2D void.

    # Step 2: Apply the Weak Morse Inequalities.
    # Morse theory states that the number of critical points of index k (c_k)
    # is greater than or equal to the k-th Betti number (b_k).
    # c_0 >= b_0 (minima)
    # c_1 >= b_1 (saddles)
    # c_2 >= b_2 (maxima)
    # Therefore, the minimal number for each type equals the Betti number.
    min_c0 = b0
    min_c1 = b1
    min_c2 = b2
    
    # Step 3: Calculate the minimal total number of critical points.
    # This is the sum of the minimal numbers for each type of critical point.
    minimal_total_points = min_c0 + min_c1 + min_c2
    
    # Step 4: Print the result as a clear equation.
    print("Based on Morse theory, the minimal number of critical points for a smooth function")
    print("on a 2-torus is determined by its Betti numbers (b_0, b_1, b_2).")
    print(f"Minimal number of minima (c_0) >= b_0 = {b0}")
    print(f"Minimal number of saddles (c_1) >= b_1 = {b1}")
    print(f"Minimal number of maxima (c_2) >= b_2 = {b2}")
    print("\nThe minimal total number of critical points is the sum:")
    print(f"{min_c0} + {min_c1} + {min_c2} = {minimal_total_points}")

solve_minimal_critical_points()
<<<4>>>