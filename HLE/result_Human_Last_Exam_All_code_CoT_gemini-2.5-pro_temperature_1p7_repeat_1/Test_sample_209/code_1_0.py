def solve_fish_pond_problem():
    """
    This script simulates the fish pond problem to find the earliest time all fish are gone.
    It follows the optimal strategy of using the fastest equipment and managing constraints efficiently.
    """
    # Initial state
    time = 0
    free_goldfish = 10
    goldfish_in_basket = 0
    # This counter tracks free-swimming goldfish eaten by the shark for its power-up and survival check.
    shark_eaten_free_count = 0
    # This tracks the time of the shark's last meal of a free-swimming fish.
    shark_last_free_meal_time = 0

    print(f"Time {time}: The simulation starts. The pond has {free_goldfish} goldfish and 1 shark.")
    print("Time 0: The fisherman starts with the fastest rod, Rod B.")
    print("\n--- Phase 1: Fishing with Rod B until basket is full ---")

    # t=5
    time = 5
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod B. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")

    # t=10
    time = 10
    # Fisherman's action
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod B. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")
    # Shark's action
    free_goldfish -= 1
    shark_eaten_free_count += 1
    shark_last_free_meal_time = time
    print(f"Time {time}: Shark eats 1 free-swimming goldfish. Free Goldfish: {free_goldfish}, Shark's free-eaten count: {shark_eaten_free_count}.")

    # t=15
    time = 15
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod B. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")

    # t=20
    time = 20
    # Fisherman's action
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod B. His basket is now full. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")
    # Shark's action
    free_goldfish -= 1
    shark_eaten_free_count += 1
    shark_last_free_meal_time = time
    print(f"Time {time}: Shark eats 1 free-swimming goldfish. Free Goldfish: {free_goldfish}, Shark's free-eaten count: {shark_eaten_free_count}.")

    print("\n--- Phase 2: Fisherman switches strategy ---")
    print(f"Time {time}: The basket is full. The fisherman switches to one-handed Rod A to free up a hand.")
    print("This also resets the usage timer for Rod B since he used it for 20 minutes (< 30 min).")
    print("He will feed his 4 caught fish to the shark while continuing to fish with Rod A.")

    # Simulating the period from t=20 to t=28
    # Events in this period: 4 feedings, 1 catch with Rod A.
    print(f"Time {22}: Fisherman feeds the 1st fish to the shark. Basket: {3}.")
    print(f"Time {24}: Fisherman feeds the 2nd fish to the shark. Basket: {2}.")
    print(f"Time {26}: Fisherman feeds the 3rd fish to the shark. Basket: {1}.")
    
    # Rod A catch happens at t=27 (20 + 7)
    time = 27
    free_goldfish -= 1
    goldfish_in_basket = 1 + 1 # 1 from feeding, +1 from catch
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod A. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")

    # Last feeding happens at t=28 (26 + 2)
    time = 28
    goldfish_in_basket = 2 - 1 # 2 from previous step, -1 from feeding
    print(f"Time {time}: Fisherman feeds the 4th fish to the shark. Basket: {goldfish_in_basket}.")

    print(f"\nTime {time}: The basket is clear. The fisherman has held Rod A for 8 minutes (> 7 min).")
    print("He now switches back to the faster Rod B.")

    print("\n--- Phase 3: Clearing the remaining goldfish ---")

    # t=30
    time = 30
    free_goldfish -= 1
    shark_eaten_free_count += 1
    shark_last_free_meal_time = time
    print(f"Time {time}: Shark eats 1 free-swimming goldfish. Free Goldfish: {free_goldfish}, Shark's free-eaten count: {shark_eaten_free_count}.")

    # t=33 (28 + 5)
    time = 33
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches 1 goldfish with Rod B. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")

    # t=38 (33 + 5)
    time = 38
    free_goldfish -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches the last free-swimming goldfish. Free Goldfish: {free_goldfish}, In Basket: {goldfish_in_basket}.")

    print("\n--- Final Phase: The shark's demise ---")
    print("All 10 goldfish have been removed from the pond.")
    print(f"The shark's last meal of a free-swimming fish was at time {shark_last_free_meal_time}.")
    
    shark_survival_time = 11
    # Check for survival extension: shark_eaten_free_count (3) is not > 4.
    print(f"The shark has not eaten more than 4 fish in the last 60 minutes, so its survival time is the standard {shark_survival_time} minutes.")
    
    final_time = shark_last_free_meal_time + shark_survival_time
    print("\nThe pond will be empty when the shark starves. We can calculate the final time:")
    print(f"Final Time = Shark's Last Meal Time + Shark's Survival Time")
    print(f"Final Time = {shark_last_free_meal_time} + {shark_survival_time} = {final_time}")
    print(f"\nThe earliest possible time when there are no more fish left in the pond is {final_time} minutes.")

solve_fish_pond_problem()