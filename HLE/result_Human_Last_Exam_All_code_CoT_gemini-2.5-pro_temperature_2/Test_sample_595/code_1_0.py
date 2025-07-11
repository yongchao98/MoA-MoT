import math

def solve():
    """
    Calculates the maximum number of grid squares a triangle with sides 18, 18, 18*sqrt(2)
    can pass through without its perimeter touching a lattice point.
    """
    N = 18

    # --- Case 1: Legs are aligned with the coordinate axes ---
    # Vertices are slightly offset from integers, e.g., A(e, e), B(18+e, e), C(e, 18+e)

    # A horizontal leg of integer length N crosses N+1 squares.
    k_leg_horz = N + 1

    # A vertical leg of integer length N crosses N+1 squares.
    k_leg_vert = N + 1

    # A diagonal (slope +/-1) leg across an N x N grid crosses 2N+1 squares.
    k_hyp_diag = 2 * N + 1

    # Overlaps (shared squares at vertices)
    # A(e,e): vertex (0,0) square, shared by the two legs.
    int_AB_AC = 1
    # B(18+e,e): vertex (18,0) square. The two squares (17,0) and (18,0) are shared.
    int_AB_BC = 2
    # C(e,18+e): vertex (0,18) square. The two squares (0,17) and (0,18) are shared.
    int_AC_BC = 2
    # No square is common to all three sides.
    int_ABC = 0

    total_k_case1 = k_leg_horz + k_leg_vert + k_hyp_diag - (int_AB_AC + int_AB_BC + int_AC_BC) + int_ABC

    print("--- Analysis ---")
    print("Case 1: Legs parallel to axes")
    print(f"Leg 1 (horizontal, length {N}) crosses {k_leg_horz} squares.")
    print(f"Leg 2 (vertical, length {N}) crosses {k_leg_vert} squares.")
    print(f"Hypotenuse (diagonal, spans {N}x{N} grid) crosses {k_hyp_diag} squares.")
    print(f"Total squares (without subtracting overlaps): {k_leg_horz + k_leg_vert + k_hyp_diag}")
    print(f"Shared squares at right-angle vertex: {int_AB_AC}")
    print(f"Shared squares at other two vertices: {int_AB_BC} and {int_AC_BC}")
    print(f"Final count for Case 1 = {k_leg_horz} + {k_leg_vert} + {k_hyp_diag} - ({int_AB_AC} + {int_AB_BC} + {int_AC_BC}) = {total_k_case1}")
    print("-" * 20)

    # --- Case 2: Rotated by 45 degrees (hypotenuse is horizontal) ---
    
    # A horizontal line of length L crosses floor(L)+1 squares.
    hyp_len = N * math.sqrt(2)
    k_hyp_horz = math.floor(hyp_len) + 1
    
    # The legs are now diagonal (slope +/- 1). The grid they span is M x M.
    # M = floor(leg_length / sqrt(2))
    M = math.floor(N / math.sqrt(2)) # M = floor(18/sqrt(2)) = floor(9*sqrt(2)) = 12
    k_leg_diag = 2 * M + 1

    # Overlaps
    # At the right-angle vertex, legs with slope +1 and -1 meet. They share 1 square.
    int_leg1_leg2 = 1
    # At the other vertices, a diagonal leg meets a horizontal line. They share 2 squares.
    int_leg_hyp_1 = 2
    int_leg_hyp_2 = 2

    total_k_case2 = k_hyp_horz + k_leg_diag + k_leg_diag - (int_leg1_leg2 + int_leg_hyp_1 + int_leg_hyp_2) + int_ABC

    print("Case 2: Rotated by 45 degrees")
    print(f"Hypotenuse (horizontal, length {hyp_len:.2f}) crosses {k_hyp_horz} squares.")
    print(f"Each leg (diagonal, spans {M}x{M} grid) crosses {k_leg_diag} squares.")
    print(f"Total squares (without subtracting overlaps): {k_hyp_horz + k_leg_diag + k_leg_diag}")
    print(f"Shared squares at right-angle vertex: {int_leg1_leg2}")
    print(f"Shared squares at other two vertices: {int_leg_hyp_1} and {int_leg_hyp_2}")
    print(f"Final count for Case 2 = {k_hyp_horz} + {k_leg_diag} + {k_leg_diag} - ({int_leg1_leg2} + {int_leg_hyp_1} + {int_leg_hyp_2}) = {total_k_case2}")
    print("-" * 20)

    final_k = max(total_k_case1, total_k_case2)
    print(f"The largest number of squares k is the maximum of the two cases.")
    print(f"k = max({total_k_case1}, {total_k_case2}) = {final_k}")

solve()