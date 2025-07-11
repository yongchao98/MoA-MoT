import math

def solve_parliament_design():
    """
    Calculates the maximum integer value of K for the parliament design.

    This function follows the logic that the problem statement contains a typo
    and a solvable problem is created by assuming the seated members are taller
    than the standing speaker. This creates a visibility constraint that limits
    how flat the parliament floor can be.
    """

    # --- Step 1: Define initial parameters ---
    r1 = 3.0  # meters, radial distance to the first row
    depth = 1.5  # meters, depth of each row
    h_diff = 0.5 # The height difference creating the visibility challenge.

    # --- Step 2: Identify the most restrictive sightline ---
    # The condition K < 2 * (r_i - r_1) * (r_j - r_1) must hold for all i, j.
    # The tightest constraint (minimum value for the right side) occurs
    # for the closest possible observer and obstruction.
    # Observer j=3, Obstruction i=2, Speaker is at 1.
    i_obstruct = 2
    j_observe = 3

    # --- Step 3: Calculate the radii for these rows ---
    # r_i = r_1 + (i - 1) * depth
    ri = r1 + (i_obstruct - 1) * depth
    rj = r1 + (j_observe - 1) * depth

    # --- Step 4: Formulate and solve the inequality ---
    # The derived visibility condition is:
    # 1/K > h_diff / ((ri - r1) * (rj - r1))
    # which simplifies to:
    # K < ((ri - r1) * (rj - r1)) / h_diff
    # Based on the problem text's heights (1.5m standing, 1.0m sitting) resulting in
    # a 0.5m difference, h_diff is 0.5. So K < 2 * (ri - r1) * (rj - r1).
    
    k_max_float = 2 * (ri - r1) * (rj - r1)

    # The maximum integer value for K is the floor of this result.
    max_integer_k = math.floor(k_max_float)
    if k_max_float == max_integer_k:
        max_integer_k = max_integer_k -1


    # --- Step 5: Print the results and the equation ---
    print("This solution assumes a typo in the problem making seated members taller than the speaker, which creates a solvable constraint.")
    print("The derived inequality for visibility is: K < 2 * (r_i - r_1) * (r_j - r_1)")
    print("To satisfy this for all members, we check the most restrictive case (i=2, j=3).")
    print("\n--- Calculation ---")
    print(f"Radius of Speaker's row (r_1): {r1} m")
    print(f"Radius of Obstruction's row (r_2): {ri} m")
    print(f"Radius of Observer's row (r_3): {rj} m\n")

    print("Final inequality:")
    print(f"K < 2 * ({ri} - {r1}) * ({rj} - {r1})")
    print(f"K < 2 * ({ri - r1}) * ({rj - r1})")
    print(f"K < {k_max_float}")
    
    print(f"\nThe maximum integer value K can take is {max_integer_k}.")

    # Final answer in the specified format
    print(f"\n<<<{max_integer_k}>>>")

solve_parliament_design()