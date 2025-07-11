def solve_horse_and_water(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter batches of water available at the start.
        m (float): The distance to the destination in kilometers.
    """
    # First, check if the journey is even possible.
    # The maximum distance is the sum of all possible stage lengths.
    max_dist = sum(100.0 / (2 * i - 1) for i in range(1, n + 1))
    if m > max_dist:
        print(f"It's impossible to travel {m} km with {n*100} liters of water.")
        print(f"The maximum possible distance is {max_dist:.2f} km.")
        return

    # 1. Find K, the number of 100L loads being transported during the final stage of the journey.
    K = 0
    dist_traveled_before_final_stage = 0.0
    # Iterate from the first stage (n loads) downwards.
    for k_current in range(n, 0, -1):
        # The distance that can be covered in a stage to consume exactly 100L of water.
        dist_of_stage = 100.0 / (2 * k_current - 1)
        
        # Check if the destination `m` falls within the current stage.
        if m <= dist_traveled_before_final_stage + dist_of_stage:
            K = k_current
            break
        else:
            # If not, we complete this stage and move to the next.
            dist_traveled_before_final_stage += dist_of_stage
            
    # 2. Build the equation strings and calculate values.
    # The summation term is the sum of distances of all stages before the final one (i.e., from n down to K+1).
    summation_parts_str = []
    for i in range(n, K, -1):
        summation_parts_str.append(f"100/(2*{i}-1)")
    
    summation_str = " + ".join(summation_parts_str) if summation_parts_str else "0"

    # 3. Print the formula and the step-by-step calculation.
    print(f"For n={n} (initial water = {n*100}L) and destination m={m}km:")
    print(f"The journey's final leg will be made while transporting K={K} loads of 100L.\n")
    
    print("The final equation for the water left is:")
    # This formula shows all the numbers that will be used in the calculation.
    final_equation_str = f"{K}*100 - ({m} - ({summation_str})) * (2*{K}-1)"
    print(final_equation_str)
    
    print("\nStep-by-step calculation:")
    
    # Step 1: Calculate the value of the summation part (distance traveled before final stage)
    sum_val = dist_traveled_before_final_stage
    print(f"= {K*100} - ({m} - {sum_val:.2f}) * (2*{K}-1)")
    
    # Step 2: Calculate the remaining distance to travel in the final stage
    remaining_dist = m - sum_val
    print(f"= {K*100} - ({remaining_dist:.2f}) * (2*{K}-1)")

    # Step 3: Calculate the total water consumed in the final stage
    water_consumed_final_stage = remaining_dist * (2 * K - 1)
    print(f"= {K*100} - {water_consumed_final_stage:.2f}")

    # Step 4: Final calculation
    final_water = K * 100 - water_consumed_final_stage
    print(f"= {final_water:.2f}\n")

    print(f"The maximum amount of water left is {final_water:.2f} liters.")
    
    # Return the final value to be captured by the system
    return final_water


if __name__ == '__main__':
    # --- User Inputs ---
    # n: The horse starts with n * 100 liters of water.
    # m: The distance to the destination in kilometers.
    n_input = 3
    m_input = 100
    # -------------------

    result = solve_horse_and_water(n_input, m_input)
    if result is not None:
        # The '<<<' and '>>>' tags are for the system to capture the final numeric answer.
        print(f"<<<{result:.2f}>>>")
