def solve_fish_pond_problem():
    """
    This script simulates the optimal strategy to clear the pond of all fish
    at the earliest possible time and calculates the final result.
    """
    
    # --- Initial State ---
    free_goldfish = 10
    goldfish_in_basket = 0
    shark_eaten_count = 0
    shark_last_meal_time = 0
    time = 0

    print("--- Phase 1: Initial Fishing with Rod B (Time 0 to 20) ---")
    print(f"Time {time}: Start. Fisherman uses Rod B. Goldfish in pond: {free_goldfish}")

    # Fisherman catches his first 4 fish, while the shark eats on its own schedule.
    # Fisherman's catches at t=5, 10, 15, 20
    # Shark's eats at t=10, 20
    
    time = 5
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")
    
    time = 10
    free_goldfish -= 1 # Fisherman's catch
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")
    free_goldfish -= 1 # Shark's meal
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats a goldfish. Free Goldfish: {free_goldfish}, Shark Eaten: {shark_eaten_count}")

    time = 15
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")

    time = 20
    free_goldfish -= 1 # Fisherman's catch
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")
    print(f"Time {time}: Basket is full. Fisherman leaves for 20 minutes to get a new one.")
    
    free_goldfish -= 1 # Shark's meal
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats a goldfish. Free Goldfish: {free_goldfish}, Shark Eaten: {shark_eaten_count}")

    print("\n--- Phase 2: Fisherman's Basket Trip (Time 20 to 40) ---")
    
    time = 30
    free_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats a goldfish. Free Goldfish: {free_goldfish}, Shark Eaten: {shark_eaten_count}")

    time = 40
    free_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats a goldfish. Free Goldfish: {free_goldfish}, Shark Eaten: {shark_eaten_count}")
    print(f"Time {time}: Fisherman returns with a new basket.")
    goldfish_in_basket = 0

    print("\n--- Phase 3: Final Fishing (Time 40 onwards) ---")
    # At t=40, there are 2 goldfish left. The fisherman continues with Rod B.
    
    time = 45
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")

    time = 50
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches the last goldfish. Free Goldfish: {free_goldfish}, Basket: {goldfish_in_basket}")
    print(f"Time {time}: All goldfish have been removed from the pond.")
    
    # The shark would have tried to eat at t=50 (10 mins after t=40), but the fisherman caught the last fish.
    # So the shark's last successful meal was at t=40.

    print("\n--- Final Calculation ---")
    
    # Check the starvation rule. At the time of its last meal (t=40), had it eaten MORE than 4 fish in the last 60 minutes?
    # It had eaten exactly 4 fish (at t=10, 20, 30, 40). 4 is not > 4.
    # Therefore, the 4-minute extension does NOT apply. The starvation time is 11 minutes.
    
    starvation_duration = 11
    print(f"The shark's last meal was at time t={shark_last_meal_time}.")
    print(f"At that time, it had eaten {shark_eaten_count} fish, which is not more than 4.")
    print(f"Therefore, its starvation time is the standard {starvation_duration} minutes.")
    
    final_time = shark_last_meal_time + starvation_duration
    
    print("\nThe earliest time when no fish are left in the pond is the time of the shark's death.")
    print("Final Time = Time of Shark's Last Meal + Shark's Starvation Duration")
    print(f"{final_time} = {shark_last_meal_time} + {starvation_duration}")
    
    print("\nThe earliest possible time is:")
    print(final_time)


solve_fish_pond_problem()
<<<51>>>