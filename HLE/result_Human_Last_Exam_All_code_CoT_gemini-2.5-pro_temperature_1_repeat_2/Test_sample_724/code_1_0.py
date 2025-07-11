import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels a distance m.

    Args:
        n (int): The number of 100-liter loads of water at the origin.
        m (float): The distance to the destination in kilometers.
    """
    if n * 100 < m:
        print("Error: Not enough water to travel the distance even with a single trip.")
        return

    initial_water = n * 100.0
    distance_covered = 0.0
    loads = n

    # These lists will store the string components for the final equation
    leg_dist_str_parts = []
    water_consumed_str_parts = []

    # --- Simulation to build the equation strings ---
    
    # Depot-laying stages (from n loads down to 2 loads)
    while loads > 1:
        consumption_rate = 2 * loads - 1
        # The maximum distance for a full stage without running out of a load
        max_leg_dist = 100.0 / consumption_rate

        if distance_covered + max_leg_dist <= m:
            # This is a full stage. We complete it.
            distance_covered += max_leg_dist
            # The distance of this leg is exactly 100 / consumption_rate
            leg_dist_str_parts.append(f"100/(2*{loads}-1)")
            # The water consumed in a full stage is exactly 100
            water_consumed_str_parts.append("100")
            loads -= 1
        else:
            # This is the final, partial leg. The journey ends here.
            # The distance of this leg is whatever is left to reach m.
            prev_dist_str = " + ".join(leg_dist_str_parts) if leg_dist_str_parts else "0"
            dist_this_leg_str = f"({m} - ({prev_dist_str}))"
            
            # Water consumed = rate * distance for this partial leg
            water_consumed_str_parts.append(f"(2*{loads}-1) * {dist_this_leg_str}")
            distance_covered = m  # Mark distance as fully covered
            break
    
    # Final one-way trip stage (if not already handled in the loop)
    if distance_covered < m:
        # This part runs if all depot stages were completed and there's still distance to cover
        # This is the final leg with only one load, so consumption rate is 1.
        prev_dist_str = " + ".join(leg_dist_str_parts) if leg_dist_str_parts else "0"
        dist_this_leg_str = f"({m} - ({prev_dist_str}))"
        
        # Water consumed = 1 * distance
        water_consumed_str_parts.append(f"1 * {dist_this_leg_str}")

    # --- Calculation of the final numerical result ---
    
    total_water_consumed = 0.0
    dist_covered_calc = 0.0
    loads_calc = n
    
    # Recalculate depot stages
    while loads_calc > 1:
        rate = 2 * loads_calc - 1
        leg_dist = 100.0 / rate
        if dist_covered_calc + leg_dist <= m:
            total_water_consumed += 100.0
            dist_covered_calc += leg_dist
            loads_calc -= 1
        else:
            rem_dist = m - dist_covered_calc
            total_water_consumed += rem_dist * rate
            dist_covered_calc = m
            break
            
    # Recalculate final one-way trip
    if dist_covered_calc < m:
        rem_dist = m - dist_covered_calc
        total_water_consumed += rem_dist * 1.0 # Consumption rate is 1

    final_water = initial_water - total_water_consumed

    # --- Print the final formatted output ---
    
    # Assemble the final equation string
    consumption_summation = " + ".join(water_consumed_str_parts)
    print(f"Initial water: {n}*100 = {initial_water:.2f} liters")
    print(f"Destination distance: {m} km")
    print("\nEquation for maximum water left at destination:")
    # The user wanted each number in the final equation. This representation shows the logic.
    equation = f"{n}*100 - ({consumption_summation})"
    print(f"{equation} = {final_water:.2f}")

    print(f"\n<<<Result>>>")
    print(f"{final_water:.2f}")


# --- User Inputs ---
# n: number of 100-liter water loads
# m: distance in kilometers
# Example values from the problem description
n_input = 3
m_input = 175

# Solve the problem with the given inputs
solve_horse_problem(n_input, m_input)

# Example 2: Short distance
# solve_horse_problem(n=4, m=20)

# Example 3: Long distance
# solve_horse_problem(n=5, m=250)