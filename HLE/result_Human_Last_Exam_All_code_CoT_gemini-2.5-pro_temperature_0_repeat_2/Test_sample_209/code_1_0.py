def solve_fishing_problem():
    """
    Calculates the earliest time the pond is empty by simulating the optimal strategy.
    The code prints a step-by-step timeline of events.
    """
    # Initial State
    time = 0
    free_goldfish = 10
    basket_goldfish = 0
    shark_last_meal_time = 0
    
    print(f"T={time:02d} min: The fisherman starts with Rod B. Pond has {free_goldfish} goldfish.")
    print("-" * 60)

    # --- Phase 1: Initial fishing with Rod B (0-15 mins) ---
    time = 5
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod B. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    
    time = 10
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod B. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    # Shark eats simultaneously
    free_goldfish -= 1
    shark_last_meal_time = 10
    print(f"T={time:02d} min: Shark eats a free-swimming goldfish. [Free: {free_goldfish}, Last Meal: T={shark_last_meal_time}]")

    time = 15
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod B. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    print("-" * 60)

    # --- Phase 2: Strategic switch to Rod A to reset Rod B timer (15-22 mins) ---
    print(f"T={time:02d} min: Fisherman switches to Rod A to reset the 30-min limit and free a hand.")
    
    time = 17 # Feed takes 2 mins
    basket_goldfish -= 1
    shark_last_meal_time = 17
    print(f"T={time:02d} min: Fisherman feeds the shark from the basket. [Basket: {basket_goldfish}, Last Meal: T={shark_last_meal_time}]")

    time = 22 # 7 mins of using Rod A have passed
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod A. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    print(f"T={time:02d} min: Rod B's timer is now reset. Fisherman switches back to Rod B.")
    print("-" * 60)

    # --- Phase 3: Back to Rod B, basket fills up (22-27 mins) ---
    time = 27 # 5 mins of using Rod B
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod B. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    # Shark eats simultaneously
    free_goldfish -= 1
    shark_last_meal_time = 27
    print(f"T={time:02d} min: Shark eats a free-swimming goldfish. [Free: {free_goldfish}, Last Meal: T={shark_last_meal_time}]")
    print(f"T={time:02d} min: Basket is now full (4 fish). Fisherman must switch to Rod A to empty it.")
    print("-" * 60)

    # --- Phase 4: Forced switch to Rod A to empty the basket (27-35 mins) ---
    time = 29 # Feed 1
    basket_goldfish -= 1
    shark_last_meal_time = 29
    print(f"T={time:02d} min: Fisherman feeds the shark. [Basket: {basket_goldfish}, Last Meal: T={shark_last_meal_time}]")
    
    time = 31 # Feed 2
    basket_goldfish -= 1
    shark_last_meal_time = 31
    print(f"T={time:02d} min: Fisherman feeds the shark. [Basket: {basket_goldfish}, Last Meal: T={shark_last_meal_time}]")

    time = 33 # Feed 3
    basket_goldfish -= 1
    shark_last_meal_time = 33
    print(f"T={time:02d} min: Fisherman feeds the shark. [Basket: {basket_goldfish}, Last Meal: T={shark_last_meal_time}]")

    time = 34 # 7 mins of Rod A passed (27 -> 34)
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod A. [Free: {free_goldfish}, Basket: {basket_goldfish}]")

    time = 35 # Feed 4
    basket_goldfish -= 1
    shark_last_meal_time = 35
    print(f"T={time:02d} min: Fisherman feeds the shark, emptying the basket. [Basket: {basket_goldfish}, Last Meal: T={shark_last_meal_time}]")
    print(f"T={time:02d} min: Rod B's timer is reset. Fisherman switches back to Rod B.")
    print("-" * 60)

    # --- Phase 5: Final fishing (35-45 mins) ---
    time = 40 # 5 mins of Rod B
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches a goldfish with Rod B. [Free: {free_goldfish}, Basket: {basket_goldfish}]")

    time = 45 # 5 more mins of Rod B
    free_goldfish -= 1
    basket_goldfish += 1
    print(f"T={time:02d} min: Fisherman catches the last free-swimming goldfish. [Free: {free_goldfish}, Basket: {basket_goldfish}]")
    print("-" * 60)

    # --- Phase 6: Shark Starvation ---
    print(f"T={time:02d} min: All goldfish have been caught. Now, we wait for the shark to starve.")
    print(f"The shark's last meal was at T={shark_last_meal_time} minutes.")
    
    base_starvation_time = 11
    # The shark ate more than 4 fish in the last 60 minutes, so its survival time is extended.
    bonus_starvation_time = 4
    total_starvation_time = base_starvation_time + bonus_starvation_time
    
    print(f"The shark's base survival time is {base_starvation_time} minutes.")
    print(f"Because it ate more than 4 fish in the last hour, it gets a {bonus_starvation_time}-minute extension.")
    
    final_time = shark_last_meal_time + total_starvation_time
    
    print("\nFinal Calculation:")
    print(f"Time of last meal + Total starvation time = Final Time")
    print(f"{shark_last_meal_time} minutes + ({base_starvation_time} + {bonus_starvation_time}) minutes = {final_time} minutes")

solve_fishing_problem()