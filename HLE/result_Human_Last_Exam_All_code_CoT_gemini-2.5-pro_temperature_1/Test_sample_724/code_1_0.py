import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter tanks of water initially available.
        m (float): The total distance to travel in km.
    """
    initial_water = n * 100

    if initial_water <= m:
        print("The horse does not have enough water to travel the distance in a single trip.")
        print("According to the problem statement (n*100 > m), this case should not occur.")
        return

    # Find k, the number of full 100-liter tanks consumed to set up depots.
    k = 0
    distance_covered_by_depots = 0.0
    summation_terms_str = []
    
    # We can set up depots as long as we have more than one tank to move.
    for i in range(1, n):
        denominator = 2 * n - 2 * i + 1
        segment_distance = 100.0 / denominator
        
        if distance_covered_by_depots + segment_distance < m:
            distance_covered_by_depots += segment_distance
            k = i
            summation_terms_str.append(f"1/{denominator}")
        else:
            break
            
    # D_k is the distance covered by setting up k depots
    D_k = distance_covered_by_depots
    
    # Water consumed to set up k depots
    water_consumed_for_depots = 100.0 * k
    
    # Remaining distance to travel
    remaining_distance = m - D_k
    
    # For the final leg, we have n-k tanks to move
    final_leg_consumption_rate = 2 * (n - k) - 1
    
    # Water consumed on the final leg
    water_consumed_final_leg = remaining_distance * final_leg_consumption_rate
    
    # Total water consumed
    total_water_consumed = water_consumed_for_depots + water_consumed_final_leg
    
    # Final amount of water left
    water_left = initial_water - total_water_consumed

    # --- Outputting the results ---
    print(f"To travel {m} km with an initial {initial_water} liters of water ({n} tanks):")
    print("-" * 50)
    
    # Express the answer in terms of summation
    print("The amount of water left is calculated by subtracting the consumed water from the initial water.")
    print("Water Consumed = (Water for Depots) + (Water for Final Leg)")
    
    if k > 0:
        sum_expr = f"100 * ({' + '.join(summation_terms_str)})"
        print(f"\nWater for Depots = 100 * k = 100 * {k} = {water_consumed_for_depots:.2f} liters")
        print(f"Distance covered by depots, D_k = {sum_expr} = {D_k:.2f} km")
        print(f"\nWater for Final Leg = (m - D_k) * (2*(n-k)-1)")
        print(f"                      = ({m} - {D_k:.2f}) * {final_leg_consumption_rate} = {water_consumed_final_leg:.2f} liters")
        
        print("\nThe final equation with all numbers is:")
        print(f"{initial_water} - ( {water_consumed_for_depots:.2f} + {water_consumed_final_leg:.2f} ) = {water_left:.2f}")

    else: # k=0, no depots are created
        print("\nNo depots are created. The entire distance is traveled with a consumption rate based on moving n tanks.")
        print(f"Water Consumed = m * (2*n - 1)")
        print(f"               = {m} * {final_leg_consumption_rate} = {total_water_consumed:.2f} liters")
        print("\nThe final equation with all numbers is:")
        print(f"{initial_water} - {total_water_consumed:.2f} = {water_left:.2f}")

    print("-" * 50)
    print(f"The maximum amount of water left is: {water_left:.2f} liters.")
    
    # The problem asks for the final answer in a specific format at the end.
    # To avoid floating point inaccuracies in the final output, we round to 2 decimal places.
    final_answer = round(water_left, 2)
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    # --- User-defined variables ---
    # n: number of 100-liter tanks
    # m: distance to destination in km
    n = 4
    m = 200

    solve_horse_problem(n, m)