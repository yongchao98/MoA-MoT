def solve_fish_pond_riddle():
    """
    Calculates the earliest possible time when there are no more fish left in the pond
    by simulating the optimal strategy for the fisherman.
    """
    # --- Define initial parameters from the problem statement ---
    initial_goldfish = 10
    basket_capacity = 4
    shark_initial_eat_interval = 10
    fisherman_rod_b_speed = 5
    basket_trip_time = 20
    shark_starvation_time = 11

    print("Calculating the earliest time for the pond to be empty...")
    print("The optimal strategy involves using the fastest tool, Rod B, and letting the shark help clear the pond.")
    print("-" * 30)

    # --- Phase 1: Fisherman fills the first basket ---
    time = 0
    fisherman_caught_phase1 = basket_capacity
    time_to_fill_basket = fisherman_caught_phase1 * fisherman_rod_b_speed
    time += time_to_fill_basket
    
    # During this time, the shark eats at regular intervals.
    shark_eaten_phase1 = time // shark_initial_eat_interval
    
    print(f"PHASE 1: Fisherman fills the first basket.")
    print(f"Equation for time taken: {fisherman_caught_phase1} goldfish * {fisherman_rod_b_speed} min/goldfish = {time_to_fill_basket} minutes.")
    print(f"In these {time_to_fill_basket} minutes, the shark eats {shark_eaten_phase1} goldfish (at t=10, t=20).")
    
    remaining_goldfish_after_phase1 = initial_goldfish - fisherman_caught_phase1 - shark_eaten_phase1
    print(f"Status at time={time} min: {remaining_goldfish_after_phase1} goldfish are left swimming freely.")
    print("-" * 30)
    
    # Keep track of the shark's last meal
    time_of_shark_last_meal = shark_eaten_phase1 * shark_initial_eat_interval

    # --- Phase 2: Fisherman travels for a new basket ---
    print(f"PHASE 2: Fisherman's {basket_trip_time}-minute trip for a new basket.")
    time_at_departure = time
    time += basket_trip_time

    shark_eats_during_trip = 0
    next_shark_meal_time = time_of_shark_last_meal + shark_initial_eat_interval
    # Calculate how many times the shark eats while the fisherman is away
    while next_shark_meal_time <= time:
        shark_eats_during_trip += 1
        time_of_shark_last_meal = next_shark_meal_time
        next_shark_meal_time += shark_initial_eat_interval
        
    print(f"During the trip (from t={time_at_departure} to t={time}), the shark eats {shark_eats_during_trip} more goldfish (at t=30, t=40).")
    remaining_goldfish_after_phase2 = remaining_goldfish_after_phase1 - shark_eats_during_trip
    total_shark_eaten = shark_eaten_phase1 + shark_eats_during_trip
    print(f"Status at time={time} min (fisherman returns): {remaining_goldfish_after_phase2} goldfish are left.")
    print(f"The shark's last meal was at t={time_of_shark_last_meal} minutes.")
    print("-" * 30)

    # --- Phase 3: Fisherman catches the remaining goldfish ---
    print(f"PHASE 3: Fisherman catches the remaining {remaining_goldfish_after_phase2} goldfish.")
    fisherman_caught_phase3 = remaining_goldfish_after_phase2
    time_to_catch_rest = fisherman_caught_phase3 * fisherman_rod_b_speed
    time_last_goldfish_caught = time + time_to_catch_rest

    print(f"Equation for time taken: {fisherman_caught_phase3} goldfish * {fisherman_rod_b_speed} min/goldfish = {time_to_catch_rest} minutes.")
    print(f"The last goldfish is caught at t={time_last_goldfish_caught} minutes.")
    print("At this point, all goldfish are gone from the pond.")
    print("-" * 30)

    # --- Phase 4: Determine the final time when the shark dies ---
    print("PHASE 4: Determine the final event.")
    print(f"The shark's last meal was at t={time_of_shark_last_meal}. There are no more fish for it to eat.")
    
    # Check if the shark's survival time is extended.
    # It ate 4 fish (at t=10, 20, 30, 40). The condition is "more than 4", so the bonus does not apply.
    survival_bonus_time = 4
    if total_shark_eaten > 4:
        shark_final_survival_time = shark_starvation_time + survival_bonus_time
        print(f"The shark ate >4 fish in the last hour, so its survival window is {shark_starvation_time} + {survival_bonus_time} = {shark_final_survival_time} minutes.")
    else:
        shark_final_survival_time = shark_starvation_time
        print(f"The shark did not eat more than 4 fish, so its survival window is the standard {shark_starvation_time} minutes.")
    
    # Calculate the final time.
    final_time = time_of_shark_last_meal + shark_final_survival_time
    print("\nFINAL CALCULATION:")
    print("The pond is empty when the shark dies. This is calculated from its last meal.")
    print(f"Final Time = (Time of Shark's Last Meal) + (Shark's Survival Time)")
    print(f"Final Time = {time_of_shark_last_meal} + {shark_final_survival_time} = {final_time}")
    
    return final_time

if __name__ == '__main__':
    final_answer = solve_fish_pond_riddle()
    # The final answer is wrapped as requested.
    # print(f"\n<<< {final_answer} >>>")