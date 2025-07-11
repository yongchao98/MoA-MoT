def solve_fish_pond_problem():
    """
    Simulates the fish pond scenario to find the earliest time all fish are gone.
    This is achieved by having the fisherman use the most efficient method available.
    """
    
    # --- Initial State ---
    time = 0
    goldfish_in_pond = 10
    
    # --- Simulation of Key Events ---
    
    print("--- Simulation Start ---")
    print("Optimal Strategy: Fisherman uses the fastest rod (Rod B) continuously.")
    print(f"Time {time}: Start. Pond has 10 goldfish. Fisherman begins fishing.")
    
    # Manually step through the determined optimal timeline
    
    time = 5
    print(f"Time {time}: Fisherman catches 1st goldfish. Pond: 9.")
    
    time = 10
    print(f"Time {time}: Fisherman catches 2nd goldfish. Pond: 8.")
    print(f"Time {time}: Shark eats 1st goldfish. Pond: 7.")
    shark_meal_times = [10]
    
    time = 15
    print(f"Time {time}: Fisherman catches 3rd goldfish. Pond: 6.")
    
    time = 20
    print(f"Time {time}: Fisherman catches 4th goldfish. Basket is now full.")
    print(f"Time {time}: Shark eats 2nd goldfish. Pond: 4.")
    shark_meal_times.append(time)
    
    print("Time 20-40: Fisherman is away for 20 minutes on a mandatory trip for a new basket.")
    
    time = 30
    print(f"Time {time}: While fisherman is away, shark eats 3rd goldfish. Pond: 3.")
    shark_meal_times.append(time)
    
    time = 40
    print(f"Time {time}: Fisherman returns. Shark eats 4th goldfish. Pond: 2.")
    shark_meal_times.append(time)
    
    time = 45
    print(f"Time {time}: Fisherman catches 5th goldfish. Pond: 1.")
    
    time = 50
    print(f"Time {time}: Fisherman catches 6th and final goldfish. Pond: 0.")
    
    time_goldfish_zero = 50
    print(f"\n--- Goldfish Depleted at Time {time_goldfish_zero} ---")
    
    # --- Final Calculation ---
    
    print("\nNow, we calculate when the shark, the last fish, is gone.")
    # Check the survival extension rule: "if the shark has eaten more than 4 fish within the last 60 minutes"
    # At t=50, the "last 60 minutes" is the window from t=-10 to t=50.
    # The shark ate 4 fish (at times 10, 20, 30, 40), all within this window.
    meals_in_last_60_mins = len(shark_meal_times)
    
    starvation_time = 11
    print(f"In the 60 minutes before the pond was empty, the shark ate {meals_in_last_60_mins} fish.")
    
    # The condition is "more than 4". Since 4 is not greater than 4, the condition is not met.
    if meals_in_last_60_mins > 4:
        starvation_time = 11 + 4
        print("This is more than 4, so the shark's survival time is extended by 4 minutes to 15.")
    else:
        print("This is not more than 4, so the standard 11-minute survival time applies.")
        
    shark_death_time = time_goldfish_zero + starvation_time
    
    print("\nThe final answer is the sum of the time when the goldfish are gone and the shark's subsequent starvation time.")
    print("Final Time = (Time Goldfish are Gone) + (Shark Starvation Time)")
    
    # Printing the numbers in the final equation as requested
    print(f"Final Time = {time_goldfish_zero} + {starvation_time} = {shark_death_time}")

# Execute the simulation and print the result.
solve_fish_pond_problem()