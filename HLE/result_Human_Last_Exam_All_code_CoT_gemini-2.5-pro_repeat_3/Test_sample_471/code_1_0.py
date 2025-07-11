def solve_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """
    # Step 1: Define the Betti numbers for the 2-torus (T^2).
    # b_0: number of connected components.
    # b_1: number of 1-dimensional "circular" holes.
    # b_2: number of 2-dimensional "voids".
    b_0 = 1
    b_1 = 2
    b_2 = 1

    print("To find the minimal number of critical points for a smooth function on a 2-torus, we use Morse theory.")
    print("The Betti numbers of the 2-torus are:")
    print(f"b_0 = {b_0} (it has 1 connected component)")
    print(f"b_1 = {b_1} (it has 2 independent loops)")
    print(f"b_2 = {b_2} (it encloses 1 void)")
    print("-" * 30)

    # Step 2 & 3: Apply the Morse inequalities and sum them.
    # The total number of critical points C must be at least the sum of the Betti numbers.
    # C = c_0 + c_1 + c_2 >= b_0 + b_1 + b_2
    print("The Morse inequalities provide a lower bound for the number of critical points:")
    print("Number of minima (c_0) >= b_0")
    print("Number of saddles (c_1) >= b_1")
    print("Number of maxima (c_2) >= b_2")
    print("\nSumming these gives the minimum total number of critical points, C_min.")
    print("C_min = b_0 + b_1 + b_2")
    print("-" * 30)

    # Step 4: Calculate the final number.
    # This bound is sharp for the torus, meaning a function exists with this exact number of critical points.
    min_critical_points = b_0 + b_1 + b_2
    
    print("The final calculation for the minimal number of critical points is:")
    # We print the final equation showing how the numbers sum up.
    print(f"{min_critical_points} = {b_0} + {b_1} + {b_2}")

solve_critical_points()