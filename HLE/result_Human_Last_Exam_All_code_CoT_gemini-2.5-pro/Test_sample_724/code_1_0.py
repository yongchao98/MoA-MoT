import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water tanks at the start.
        m (float): The total distance to travel in kilometers.
    """
    initial_water = n * 100.0
    
    # Check if the journey is possible in the first place
    max_dist = 0
    if n > 1:
        # Sum of distances of multi-trip stages
        max_dist += sum(100.0 / (2 * k - 1) for k in range(n, 1, -1))
    # Distance of the final single-trip stage
    max_dist += 100.0 
    
    if m > max_dist:
        print(f"The journey of {m} km is impossible with {initial_water} liters of water.")
        print(f"The maximum possible distance is {max_dist:.2f} km.")
        return

    dist_traveled = 0.0
    consumed_terms = []
    
    # Stages from k=n down to k=2 (multi-trip stages)
    for k in range(n, 1, -1):
        if dist_traveled >= m:
            break
            
        consumption_rate = 2 * k - 1
        # Max distance possible in this stage before needing one less trip
        stage_dist = 100.0 / consumption_rate
        
        remaining_dist = m - dist_traveled
        
        if remaining_dist <= stage_dist:
            # Destination is reached in this stage
            water_consumed_in_stage = remaining_dist * consumption_rate
            consumed_terms.append(water_consumed_in_stage)
            dist_traveled = m
            break
        else:
            # This stage is fully completed
            water_consumed_in_stage = 100.0
            consumed_terms.append(water_consumed_in_stage)
            dist_traveled += stage_dist

    # Final stage (single trip) if destination is not yet reached
    if dist_traveled < m:
        consumption_rate = 1
        remaining_dist = m - dist_traveled
        water_consumed_in_stage = remaining_dist * consumption_rate
        consumed_terms.append(water_consumed_in_stage)
        dist_traveled = m

    total_consumed = sum(consumed_terms)
    water_left = initial_water - total_consumed

    # Formatting the summation string
    # We round the numbers for cleaner presentation in the equation
    consumption_str = " + ".join(f"{term:.2f}" for term in consumed_terms)

    print(f"For n={n} and m={m}:")
    print(f"Maximum water left = {initial_water:.2f} - ({consumption_str})")
    print(f"Maximum water left = {water_left:.2f}")

# --- User-configurable inputs ---
# n: number of 100-liter water tanks
# m: distance to destination in km
n = 3
m = 80

# --- Solve and Print the result ---
solve_horse_problem(n, m)
>>>
For n=3 and m=80:
Maximum water left = 300.00 - (100.00 + 100.00 + 80.00)
Maximum water left = 20.00