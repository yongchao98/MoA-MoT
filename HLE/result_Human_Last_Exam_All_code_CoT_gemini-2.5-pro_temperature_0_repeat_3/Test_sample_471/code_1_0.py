def solve_critical_points_torus():
    """
    Calculates the minimal number of critical points for a smooth Morse function
    on a 2-torus using Morse theory.
    """
    # Step 1: Define the Betti numbers for the 2-torus (T^2).
    # b0 is the number of connected components.
    # b1 is the number of 1-dimensional "circular" holes.
    # b2 is the number of 2-dimensional "voids".
    b0 = 1
    b1 = 2
    b2 = 1

    print(f"The topology of the 2-torus is described by its Betti numbers:")
    print(f"b0 = {b0}, b1 = {b1}, b2 = {b2}\n")

    # Step 2: Use the Morse inequalities (c_k >= b_k) to find the minimal
    # number of critical points of each type for a Morse function.
    # c0: number of minima (index 0)
    # c1: number of saddle points (index 1)
    # c2: number of maxima (index 2)
    c0_min = b0
    c1_min = b1
    c2_min = b2

    print("From the Morse inequalities, the minimal number of critical points of each type are:")
    print(f"Minimal minima (c0) = {c0_min}")
    print(f"Minimal saddles (c1) = {c1_min}")
    print(f"Minimal maxima (c2) = {c2_min}\n")

    # Step 3: Calculate the total minimal number of critical points.
    total_critical_points = c0_min + c1_min + c2_min

    print("The minimal total number of critical points is the sum of these minima:")
    print(f"Equation: {c0_min} + {c1_min} + {c2_min} = {total_critical_points}\n")

    # Step 4: Verify with the Euler-Poincaré formula.
    # The formula states: c0 - c1 + c2 = chi(T^2), where chi is the Euler characteristic.
    # chi(T^2) = b0 - b1 + b2
    chi = b0 - b1 + b2
    
    print("This result is consistent with the Euler-Poincaré formula (c0 - c1 + c2 = chi):")
    print(f"Euler characteristic (chi) = {b0} - {b1} + {b2} = {chi}")
    print(f"Check: {c0_min} - {c1_min} + {c2_min} = {c0_min - c1_min + c2_min}")
    print(f"Since {c0_min - c1_min + c2_min} equals {chi}, the numbers are consistent.\n")

    print(f"Thus, the minimal number of critical points for a smooth Morse function on a 2-torus is {total_critical_points}.")

solve_critical_points_torus()