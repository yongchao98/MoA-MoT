def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function on a 2-torus
    using Morse theory.
    """
    # Step 1: Define the Betti numbers for the 2-torus (T^2).
    # b0 is the number of connected components. A torus is one connected piece.
    b0 = 1
    # b1 is the number of 1-dimensional "holes". A torus has two independent loops
    # (one around its "tube" and one through its "hole").
    b1 = 2
    # b2 is the number of 2-dimensional "voids". A torus encloses a single void.
    b2 = 1

    print("The Betti numbers for the 2-torus (T^2) are:")
    print(f"b_0 (connected components) = {b0}")
    print(f"b_1 (handles/tunnels) = {b1}")
    print(f"b_2 (voids) = {b2}\n")

    # Step 2: Apply the Morse inequalities.
    # The number of critical points of index k, c_k, must satisfy c_k >= b_k.
    # This means there must be at least:
    # - c_0 >= b_0 = 1 minimum (index 0)
    # - c_1 >= b_1 = 2 saddle points (index 1)
    # - c_2 >= b_2 = 1 maximum (index 2)
    # The minimal number of critical points for a Morse function is the sum of the Betti numbers.
    # This result holds for any smooth function in general.
    
    # These are the minimal number of critical points of each type (index).
    min_c0 = b0
    min_c1 = b1
    min_c2 = b2
    
    # Step 3: Calculate the total minimal number of critical points.
    min_total_critical_points = min_c0 + min_c1 + min_c2

    print("According to Morse theory, a smooth function on the torus must have at least:")
    print(f"- {min_c0} minimum (index 0)")
    print(f"- {min_c1} saddle points (index 1)")
    print(f"- {min_c2} maximum (index 2)\n")
    
    print("The minimal total number of critical points is the sum of these values:")
    # The final output prints the equation for the sum
    print(f"{min_c0} + {min_c1} + {min_c2} = {min_total_critical_points}")

solve_minimal_critical_points()