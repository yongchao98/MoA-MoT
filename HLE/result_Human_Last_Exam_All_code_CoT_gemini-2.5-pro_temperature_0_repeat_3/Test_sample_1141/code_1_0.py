import math

def solve_cutting_problem():
    """
    Analyzes and solves the given 3D cutting stock problem based on the user's formulation.
    """

    # --- Analysis of the Formulation ---
    # The formulation has a critical issue with the non-overlapping constraints involving T1 cubes.
    # Specifically, the constraint for a T1 cube and a B2 ball is:
    # min(|x_t - x_b|, |y_t - y_b|, |z_t - z_b|) >= 5
    #
    # A B2 ball's center must have z_b = 4.
    # A T1 cube's center must have z_t in the range [1, 7].
    # The constraint on the z-axis requires |z_t - 4| >= 5.
    # - If z_t > 4, then z_t must be >= 9. This is outside the allowed range [1, 7].
    # - If z_t <= 4, then z_t must be <= -1. This is also outside the allowed range.
    #
    # Conclusion: According to your formulation, it is impossible to place even one T1 cube
    # and one B2 ball in the same billet. Therefore, we must find the best solution by
    # considering two separate scenarios and choosing the one with the higher value.

    # --- Scenario 1: Pack only T1 cubes and B1 balls ---
    # The T1-T1 non-overlapping constraint is: min(|x_i - x_j|, |y_i - y_j|, |z_i - z_j|) >= 2.
    # This means that for any two T1 cubes, their centers must be separated by at least 2 units
    # in *every* dimension.
    # Let's define the possible "slots" for center coordinates based on this:
    # X_slots: {1, 3, 5, ..., 31} (16 slots)
    # Y_slots: {1, 3, 5, ..., 21} (11 slots)
    # Z_slots: {1, 3, 5, 7}      (4 slots)
    # To satisfy the constraint, each T1 must have a unique x, y, and z slot. The maximum
    # number of such items is limited by the smallest set of slots, which is Z_slots.
    # Max T1s = min(16, 11, 4) = 4.
    # The T1-B1 constraint is the same, so if we place 4 T1s, we use all 4 z-slots,
    # leaving no room for B1s.
    value_scenario1 = 4 * 5 + 0 * 1  # 4 T1s, 0 B1s

    # --- Scenario 2: Pack only B2 balls and B1 balls ---
    # The value of a B2 (150) is much higher than a B1 (1), so we prioritize packing B2s.
    # B2-B2 constraint: (x_i-x_j)^2 + (y_i-y_j)^2 >= 64 (since z is always 4).
    # This means the 2D distance between centers in the XY plane must be at least 8.
    # We can place 8 B2s in a 4x2 grid pattern:
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (4, 12, 4), (12, 12, 4), (20, 12, 4), (28, 12, 4)
    ]
    num_b2 = len(b2_centers)

    # Now, we perform a greedy search to fill the remaining space with B1 balls.
    b1_centers = []
    # Iterate through all possible B1 center locations.
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                new_b1_center = (x, y, z)
                is_valid = True

                # Check non-overlap with B2s: dist_sq >= 25
                for b2_c in b2_centers:
                    dist_sq = (new_b1_center[0] - b2_c[0])**2 + \
                              (new_b1_center[1] - b2_c[1])**2 + \
                              (new_b1_center[2] - b2_c[2])**2
                    if dist_sq < 25:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Check non-overlap with already placed B1s: dist_sq >= 4
                for b1_c in b1_centers:
                    dist_sq = (new_b1_center[0] - b1_c[0])**2 + \
                              (new_b1_center[1] - b1_c[1])**2 + \
                              (new_b1_center[2] - b1_c[2])**2
                    if dist_sq < 4:
                        is_valid = False
                        break
                
                if is_valid:
                    b1_centers.append(new_b1_center)
    
    num_b1 = len(b1_centers)
    value_scenario2 = num_b2 * 150 + num_b1 * 1

    # --- Conclusion ---
    print("Analysis of the Problem Formulation:")
    print("The provided formulation is unusual. The non-overlapping constraints for T1 cubes are overly restrictive and make it impossible to combine T1 cubes with B2 balls.")
    print("\nSolving the problem as formulated requires choosing the best of two scenarios:\n")
    print(f"Scenario A (T1s and B1s): The maximum value is ${value_scenario1}.")
    print(f"Scenario B (B2s and B1s): The maximum value is found by packing {num_b2} B2 balls and then filling the remaining space with B1 balls.")
    
    print("\nFinal Calculation:")
    print("The highest value is achieved in Scenario B.")
    
    final_value = max(value_scenario1, value_scenario2)
    
    # Output the final equation with each number
    print(f"{num_b2} * 150 + {num_b1} * 1 = {final_value}")

if __name__ == '__main__':
    solve_cutting_problem()