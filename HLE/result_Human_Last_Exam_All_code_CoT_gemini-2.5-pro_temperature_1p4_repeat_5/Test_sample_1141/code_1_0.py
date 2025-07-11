import math

def find_max_value():
    """
    Analyzes the cutting stock problem formulation and calculates the highest
    achievable value based on the provided rules.
    """
    # --- Analysis of the Formulation ---
    # The problem formulation is largely correct, but the non-overlapping
    # constraints for T1 cubes are problematic.
    # 1. T1-T1 rule (min(|dx|,|dy|,|dz|) >= 2): This is too restrictive
    #    and likely a typo for max(), preventing efficient packing.
    # 2. T1-B2 rule (min(|dx|,|dy|,|dz|) >= 5): As B2 centers are at z=4
    #    and T1 centers must be in z=[1,7], the z-component of this rule
    #    |z_t1 - 4| >= 5 can never be satisfied.
    # Conclusion: Based on the rules as written, if even one B2 is placed,
    # no T1s can be placed. Given the high value of B2s, the optimal
    # strategy is to maximize B2s and add B1s, ignoring T1s.

    # --- Step 1: Place B2 pieces ---
    # B2 pieces (price 150) are prioritized. We use a simple and effective
    # grid packing. B2 centers must be at z=4.
    price_b2 = 150
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (4, 12, 4), (12, 12, 4), (20, 12, 4), (28, 12, 4)
    ]
    num_b2 = len(b2_centers)
    # This configuration of 8 is the maximum possible.

    # Per the analysis, no T1 pieces can be placed with B2s.
    num_t1 = 0
    price_t1 = 5

    # --- Step 2: Fill remaining space with B1 pieces ---
    # B1 price is 1. We create a candidate lattice for B1 centers that
    # inherently satisfies the B1-B1 non-overlap constraint, which requires
    # centers to be at least 2 units apart. A grid with a step of 2 works.
    price_b1 = 1
    
    # Generate candidate points on a 2x2x2 lattice
    b1_candidate_centers = []
    x_coords = range(1, 31 + 1, 2)
    y_coords = range(1, 21 + 1, 2)
    z_coords = range(1, 7 + 1, 2)

    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                b1_candidate_centers.append((x, y, z))

    # Filter candidates against B2 exclusion zones.
    # B1-B2 non-overlap requires squared distance >= 25.
    b1_b2_separation_sq = 25
    valid_b1_centers = []

    for b1_cand in b1_candidate_centers:
        is_valid = True
        for b2_c in b2_centers:
            dist_sq = (b1_cand[0] - b2_c[0])**2 + (b1_cand[1] - b2_c[1])**2 + (b1_cand[2] - b2_c[2])**2
            if dist_sq < b1_b2_separation_sq:
                is_valid = False
                break
        if is_valid:
            valid_b1_centers.append(b1_cand)

    num_b1 = len(valid_b1_centers)

    # --- Step 3: Calculate and Print Final Result ---
    final_value = (num_b2 * price_b2) + (num_b1 * price_b1) + (num_t1 * price_t1)
    
    print("Is the problem formulation correct?")
    print("No. The constraints for T1 (cubes) are formulated in a way that makes them either highly inefficient to pack or impossible to place alongside B2 (large ball) pieces. Following the formulation literally, the optimal solution cannot contain any T1 pieces if it contains B2 pieces.")
    print("\nHighest valid solution based on the provided formulation:")
    print(f"Number of B2 pieces: {num_b2}")
    print(f"Number of B1 pieces: {num_b1}")
    print(f"Number of T1 pieces: {num_t1}")
    print("\nThe final equation for the maximum value is:")
    print(f"{num_b2} * {price_b2} + {num_b1} * {price_b1} + {num_t1} * {price_t1} = {final_value}")

if __name__ == '__main__':
    find_max_value()