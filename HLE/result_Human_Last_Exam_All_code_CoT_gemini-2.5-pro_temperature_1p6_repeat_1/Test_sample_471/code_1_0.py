def solve_critical_points_on_torus():
    """
    Calculates the minimal number of critical points for a smooth function on a 2-torus
    using Morse theory.
    """
    print("Step 1: Understanding the topology of the 2-torus (T^2) using Betti numbers.")
    # Betti numbers for the 2-torus T^2 = S^1 x S^1
    # b_0: number of connected components
    b0_torus = 1
    # b_1: number of 1D or "circular" holes
    b1_torus = 2
    # b_2: number of 2D or "cavity" holes
    b2_torus = 1
    print(f"The Betti numbers for the 2-torus are:")
    print(f"  b_0 = {b0_torus} (it's one connected piece)")
    print(f"  b_1 = {b1_torus} (it has two independent loops, like latitude and longitude)")
    print(f"  b_2 = {b2_torus} (it encloses one cavity)")
    print("-" * 20)

    print("Step 2: Applying the Morse Inequalities.")
    print("Morse theory states that for a smooth Morse function on a manifold, the number")
    print("of critical points of index k (c_k) must be at least the k-th Betti number (b_k).")
    print("So, c_k >= b_k.")
    print("\nFor the 2-torus, the critical points are classified by index:")
    print("  c_0: number of local minima")
    print("  c_1: number of saddle points")
    print("  c_2: number of local maxima")
    print("")

    # Applying the inequalities
    min_c0 = b0_torus
    min_c1 = b1_torus
    min_c2 = b2_torus
    print(f"The minimal number of minima (c_0) must be at least b_0: c_0 >= {b0_torus}")
    print(f"The minimal number of saddles (c_1) must be at least b_1: c_1 >= {b1_torus}")
    print(f"The minimal number of maxima (c_2) must be at least b_2: c_2 >= {b2_torus}")
    print("-" * 20)
    
    print("Step 3: Calculating the total minimal number of critical points.")
    print("The minimal total is the sum of the minimal values for each index.")
    
    # Calculate the total
    total_min_critical_points = min_c0 + min_c1 + min_c2
    
    # Print the final equation as requested
    print("\nFinal Calculation:")
    print(f"Minimal total = (min c_0) + (min c_1) + (min c_2)")
    print(f"              = {min_c0} + {min_c1} + {min_c2} = {total_min_critical_points}")
    print("\nThis theoretical minimum of 4 can be achieved. For example, the function")
    print("f(x, y) = cos(x) + cos(y) on the torus (defined by x, y in [0, 2pi])")
    print("has exactly 1 minimum, 2 saddle points, and 1 maximum.")

solve_critical_points_on_torus()
