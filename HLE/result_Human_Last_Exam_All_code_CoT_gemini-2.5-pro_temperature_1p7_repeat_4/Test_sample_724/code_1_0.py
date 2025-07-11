def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.
    It prints the step-by-step equation and the final numerical result.

    Args:
        n (int): The number of 100-liter water loads at the origin.
        m (float or int): The total distance to the destination in km.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(m, (int, float)) or m < 0:
        print("Error: m must be a non-negative number.")
        return
    if n * 100 < m:
        print("Error: The trip is impossible. The total water is less than the distance.")
        return
        
    # Handle the trivial case of zero distance
    if m == 0:
        initial_water = n * 100.0
        print(f"Since the distance is 0, no water is consumed.")
        print(f"Maximum water left = {initial_water:.4f}")
        print(f"Equation: {n}*100 - 0")
        print(f"<<<{initial_water:.4f}>>>")
        return

    # Step 1: Find which stage the journey ends in.
    # A stage is defined by the number of 100L loads ('k') being transported.
    distance_covered_before_stage = 0.0
    final_stage_k = 0  # k represents the number of loads being moved in a stage.

    # Iterate from n loads down to 2 loads (multi-trip stages)
    for k in range(n, 1, -1):
        consumption_rate = 2 * k - 1
        max_dist_in_stage = 100.0 / consumption_rate
        
        if distance_covered_before_stage + max_dist_in_stage >= m:
            final_stage_k = k
            break
        
        distance_covered_before_stage += max_dist_in_stage

    # If the loop completes, the destination is in the final, single-trip stage (k=1)
    if final_stage_k == 0 and n >= 1:
        final_stage_k = 1

    # Step 2: Construct the equation string part by part.
    
    # Part A: Sum of water consumed in the fully completed stages.
    num_full_stages = n - final_stage_k
    full_stages_consumption_str_parts = []
    if num_full_stages > 0:
        # Each full stage consumes exactly 100L.
        consumption_str = " + ".join(["100"] * num_full_stages)
        if num_full_stages > 1:
             full_stages_consumption_str_parts.append(f"({consumption_str})")
        else:
             full_stages_consumption_str_parts.append(consumption_str)

    # Part B: Sum of distances of the fully completed stages.
    full_stages_dist_str_parts = []
    for i in range(n, final_stage_k, -1):
        full_stages_dist_str_parts.append(f"100/(2*{i}-1)")
    
    if not full_stages_dist_str_parts:
        full_stages_dist_sum_str = "0"
    else:
        dist_sum_str = " + ".join(full_stages_dist_str_parts)
        # Add parentheses for clarity if there's more than one term
        if len(full_stages_dist_str_parts) > 1:
            full_stages_dist_sum_str = f"({dist_sum_str})"
        else:
            full_stages_dist_sum_str = dist_sum_str

    # Part C: Water consumed in the final (potentially partial) stage.
    # Formula: (m - distance_of_full_stages) * consumption_rate_of_final_stage
    rem_dist_str = f"({m} - {full_stages_dist_sum_str})"
    final_stage_consumption_str = f"{rem_dist_str} * (2*{final_stage_k}-1)"

    # Combine all consumption parts into a single string for the final equation.
    total_consumption_parts = full_stages_consumption_str_parts + [final_stage_consumption_str]
    total_consumption_str = " + ".join(total_consumption_parts)
    
    # Final equation for the water left at the destination.
    initial_water_str = f"{n}*100"
    equation = f"{initial_water_str} - ({total_consumption_str})"
    
    # Step 3: Calculate the final numerical result by evaluating the built equation string.
    # eval() is safe here because we have constructed the string ourselves from numbers and basic operators.
    water_left = eval(equation)

    # Print the results in the required format.
    print(f"For n={n} loads ({n*100} L) and a distance m={m} km:")
    print("\nThe maximum amount of water left is calculated by subtracting the total water consumed from the initial amount.")
    print("The total water consumed is the sum of water consumed in each stage of the journey.")
    print("\nThe calculation can be expressed as:")
    print(equation)
    print(f"\nFinal Answer: {water_left:.4f} liters of water are left.")
    
    # Output the final answer in the specified format
    print(f"<<<{water_left:.4f}>>>")


if __name__ == '__main__':
    # --- Set your input values here ---
    n = 4       # The horse starts with n*100 liters of water
    m = 100     # The destination is m kilometers away
    # ------------------------------------
    
    solve_horse_problem(n, m)
