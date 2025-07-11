def solve_horse_water_problem():
    """
    Calculates the maximum amount of water left after a horse journey,
    using an optimal strategy of creating water depots.
    """
    # --- User-configurable values ---
    # n: represents n*100 liters of initial water
    n = 3
    # m: the distance to the destination in kilometers
    m = 100
    # ------------------------------------

    # Per the problem, we assume n*100 > m, so a simple trip is always possible.
    # The optimal strategy, however, yields more water at the end.

    # Step 1: Find 'p', the number of fully completed transport segments.
    p = 0
    distance_of_p_segments = 0.0
    for k in range(1, n):
        # Consumption rate for segment k, determined by the number of round trips needed.
        consumption_rate = 2 * (n - (k - 1)) - 1
        
        # Length of a segment is the distance covered to consume exactly 100L of water.
        segment_length = 100.0 / consumption_rate

        if distance_of_p_segments + segment_length <= m:
            distance_of_p_segments += segment_length
            p = k
        else:
            # The destination 'm' is reached within the current segment.
            break

    # Step 2: Build the equation string and calculate the final result.
    
    # Build the string for the summation part of the formula.
    # This represents the total distance covered by the 'p' completed segments.
    sum_terms = []
    for k in range(1, p + 1):
        rate = 2 * n - 2 * k + 1
        sum_terms.append(f"100/({rate})")

    # If p=0, no segments were completed, so the sum is 0.
    sum_str = " + ".join(sum_terms) if p > 0 else "0"

    # Water available at the start of the final, partial segment.
    water_after_p_segments = (n - p) * 100.0
    
    # Consumption rate for this final leg of the journey.
    rate_of_final_leg = 2 * (n - p) - 1
    
    # Final calculation for the remaining water.
    final_water_left = water_after_p_segments - (m - distance_of_p_segments) * rate_of_final_leg
    
    # Output the final equation with all numbers substituted, showing the calculation.
    print("Equation for maximum water left:")
    print(f"({n}-{p})*100 - ({m} - ({sum_str})) * {rate_of_final_leg} = {final_water_left:.4f}")

solve_horse_water_problem()