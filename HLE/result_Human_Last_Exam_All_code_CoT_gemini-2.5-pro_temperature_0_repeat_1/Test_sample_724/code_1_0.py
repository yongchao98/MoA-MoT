def solve_horse_water_problem():
    """
    Calculates the maximum amount of water left after a horse travels a specified
    distance with a given initial amount of water, based on the optimal depot strategy.
    """
    # --- Configuration ---
    # n: the number of 100-liter tanks the horse starts with.
    # m: the total distance in kilometers to the destination.
    n = 5  # Example: 5 * 100 = 500 liters
    m = 180 # Example: 180 km

    # --- Initialization ---
    initial_water = n * 100
    distance_left = float(m)
    
    # This list will store the amount of water consumed in each phase of the journey.
    consumed_parts = []

    print(f"Calculating water left for a journey of {m} km with an initial {initial_water} liters of water.\n")

    # --- Phase 1: Multi-trip legs ---
    # Loop from 'n' tanks down to 2. Each iteration represents a phase where
    # the horse moves the entire remaining water supply forward.
    for current_tanks in range(n, 1, -1):
        if distance_left <= 0:
            break

        # Cost in liters per km for this phase.
        # (current_tanks) forward trips + (current_tanks - 1) return trips.
        cost_per_km = 2 * current_tanks - 1
        
        # The maximum distance that can be covered in this phase by consuming exactly 100L of water.
        distance_for_full_phase = 100.0 / cost_per_km

        if distance_left > distance_for_full_phase:
            # If the remaining distance is greater than what this phase covers,
            # we complete the full phase, consuming exactly 100L.
            water_consumed = 100.0
            consumed_parts.append(f"{water_consumed:.0f}")
            distance_left -= distance_for_full_phase
        else:
            # If the destination is reached within this phase, calculate consumption
            # for only the remaining distance and then stop.
            water_consumed = distance_left * cost_per_km
            consumed_parts.append(f"{water_consumed:.2f}")
            distance_left = 0
            break
            
    # --- Phase 2: Final one-way leg ---
    # If there is still distance to cover, the horse is on its last 100L tank
    # and travels one-way.
    if distance_left > 0:
        water_consumed = distance_left
        consumed_parts.append(f"{water_consumed:.2f}")

    # --- Final Calculation and Output ---
    total_consumed = sum(float(part) for part in consumed_parts)
    final_water_amount = initial_water - total_consumed

    # Build the equation string for clear representation
    equation_str = f"{initial_water}"
    if consumed_parts:
        equation_str += " - " + " - ".join(consumed_parts)

    print("The maximum amount of water left is calculated by the equation:")
    print(f"{equation_str} = {final_water_amount:.2f} liters")

    if final_water_amount < 0:
        print("\nNote: The result is negative, which means the horse cannot reach the destination.")
    
    # Print the final answer in the specified format
    print(f"\n<<<{final_water_amount:.2f}>>>")

# Execute the function
solve_horse_water_problem()