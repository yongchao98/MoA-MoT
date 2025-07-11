def solve_fish_problem():
    """
    This script calculates the earliest time all fish are gone from the pond
    by simulating the optimal strategy for the fisherman.
    """
    # --- Initial State ---
    time = 0
    pond_goldfish = 10
    basket_goldfish = 0
    shark_eaten_total = 0
    shark_last_meal_time = 0
    
    print("--- Fish Problem Solution ---")
    print(f"Initial state at T=0 min: {pond_goldfish} goldfish in pond, {basket_goldfish} in basket.\n")

    # --- Phase 1: Fisherman uses Rod B (catches 1 fish/5 min) ---
    # The fisherman uses the faster Rod B. The shark eats on its 10-minute schedule.
    # This phase continues until the basket is full (4 fish).
    print("--- Phase 1: Fisherman uses Rod B to fill the basket ---")
    
    # Event at T=5
    time = 5
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time} min: Fisherman catches a goldfish with Rod B. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")

    # Events at T=10
    time = 10
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time} min: Fisherman catches a goldfish with Rod B. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")
    pond_goldfish -= 1
    shark_eaten_total += 1
    shark_last_meal_time = time
    print(f"T={time} min: Shark eats a free-swimming goldfish. Pond: {pond_goldfish}, Shark has eaten: {shark_eaten_total}.")

    # Event at T=15
    time = 15
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time} min: Fisherman catches a goldfish with Rod B. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")

    # Events at T=20
    time = 20
    pond_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time} min: Fisherman catches a goldfish with Rod B. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")
    print("          Basket is now full. The fisherman must empty it.")
    pond_goldfish -= 1
    shark_eaten_total += 1
    shark_last_meal_time = time
    print(f"T={time} min: Shark eats a free-swimming goldfish. Pond: {pond_goldfish}, Shark has eaten: {shark_eaten_total}.\n")

    # --- Phase 2: Fisherman feeds the shark to empty the basket ---
    # This is faster (8 mins) than getting a new basket (20 mins).
    # The fisherman switches to Rod A to free a hand for feeding.
    print("--- Phase 2: Fisherman switches to Rod A and feeds the shark ---")
    
    rod_a_start_time = time
    
    # Fisherman feeds the 4 fish from the basket (2 mins each).
    # During this, he also catches one fish with Rod A at T=27.
    print(f"T={time} min: Starts feeding the 4 caught fish to the shark.")
    time_of_rod_a_catch = rod_a_start_time + 7
    
    for i in range(1, 5):
        time += 2
        basket_goldfish -= 1
        shark_eaten_total += 1
        shark_last_meal_time = time
        print(f"T={time} min: Feeds shark fish #{i}. Basket: {basket_goldfish}, Shark has eaten: {shark_eaten_total}.")
        
        if time > time_of_rod_a_catch and basket_goldfish == 1: # Check if the catch happened
             pond_goldfish -= 1
             basket_goldfish += 1
             print(f"T={time_of_rod_a_catch} min: Catches a goldfish with Rod A. Pond: {pond_goldfish}, Basket: {basket_goldfish}.")

    print(f"\nAt T={time} min, the shark has eaten {shark_eaten_total} fish and becomes stronger (eats every 2 mins).\n")
    
    # --- Phase 3: Empowered shark clears the pond ---
    print("--- Phase 3: The empowered shark eats the remaining goldfish ---")
    time_goldfish_gone = time
    while pond_goldfish > 0:
        time_goldfish_gone += 2
        pond_goldfish -= 1
        shark_eaten_total += 1
        shark_last_meal_time = time_goldfish_gone
        print(f"T={time_goldfish_gone} min: Shark eats a goldfish. Pond: {pond_goldfish}.")

    print(f"\nAt T={time_goldfish_gone} min, all goldfish are gone from the pond.\n")

    # --- Phase 4: Shark starvation ---
    print("--- Phase 4: Calculate when the shark starves ---")
    base_starvation_time = 11
    bonus_starvation_time = 4
    # The shark ate > 4 fish in the last 60 mins, so it gets the bonus.
    total_starvation_time = base_starvation_time + bonus_starvation_time
    
    print(f"The shark's last meal was at T={shark_last_meal_time} min.")
    print(f"Its starvation time is {base_starvation_time} (base) + {bonus_starvation_time} (bonus) = {total_starvation_time} minutes.")
    
    final_time = shark_last_meal_time + total_starvation_time
    
    print("\nThe final time is when the shark starves after all goldfish are gone.")
    print("Final Time = (Time when last goldfish is gone) + (Shark starvation time)")
    print(f"Final Time = {shark_last_meal_time} + {total_starvation_time} = {final_time}")
    
    print(f"\n<<<The earliest possible time when there are no more fish left in the pond is {final_time} minutes.>>>")

solve_fish_problem()