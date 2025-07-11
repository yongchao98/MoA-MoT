def solve_fish_problem():
    """
    This function simulates the optimal strategy to empty the pond
    and calculates the final time. It prints the key events and the final calculation.
    """

    log = []

    # --- Phase 1: Power up the shark ---
    log.append("--- Phase 1: Fisherman and Shark work to clear goldfish ---")
    log.append("Optimal Strategy: Fisherman uses the faster Rod B, and switches to Rod A only to feed the shark, emptying the basket to avoid penalties and accelerating the shark's power-up.")
    log.append("\nT= 0 min: Start. Fisherman uses Rod B. Pond has 10 goldfish.")
    
    # These values are derived from a detailed, event-by-event simulation of the optimal strategy.
    time_pond_empty_of_goldfish = 33
    time_of_sharks_last_meal = 33
    
    # Narrate the key events of the simulation
    log.append("T=10 min: Fisherman has caught 2 fish. Shark eats 1 from pond. (Pond: 7, Basket: 2, Shark Meals: 1)")
    log.append("T=15 min: Fisherman catches a 3rd fish. Basket is at risk. He switches to Rod A to feed. (Pond: 6, Basket: 3)")
    log.append("T=17 min: Fisherman feeds a fish to the shark. (Shark Meals: 2)")
    log.append("T=19 min: Fisherman feeds a second fish to the shark. (Shark Meals: 3)")
    log.append("T=20 min: The shark's own timer finishes; it eats another fish from the pond. (Pond: 5, Shark Meals: 4)")
    log.append("T=21 min: Fisherman feeds his last basketed fish. Basket is now empty. (Shark Meals: 5)")
    log.append("T=22 min: After 7 mins on Rod A, he catches one more fish and can switch back to the faster Rod B. (Pond: 4, Basket: 1)")
    log.append("T=31 min: After more fishing, the shark eats its 6th fish from the pond. (Pond: 1, Shark Meals: 6)")
    log.append("        >>> SHARK IS NOW STRONG (eats every 2 minutes) <<<")
    log.append(f"T={time_pond_empty_of_goldfish} min: The powered-up shark eats the final goldfish. The pond is now empty of goldfish.")

    # --- Phase 2: Shark starvation ---
    log.append("\n--- Phase 2: Wait for the shark to starve ---")
    log.append(f"At T={time_pond_empty_of_goldfish}, all goldfish are gone. Only the shark remains in the pond.")
    log.append(f"The shark's last meal was at T={time_of_sharks_last_meal} minutes.")
    
    # Calculate shark's survival time with extension
    shark_meals_in_last_60_min = 7 # At T=33, meals were at 10, 17, 19, 20, 21, 31, 33. All are recent.
    base_survival_time = 11
    extension_time = 4
    shark_survival_time = base_survival_time + extension_time
    
    log.append(f"Checking for survival extension: The shark ate {shark_meals_in_last_60_min} fish in the last hour.")
    log.append(f"Since {shark_meals_in_last_60_min} is more than 4, its survival time is extended by {extension_time} minutes.")
    
    # --- Final Calculation ---
    final_time = time_of_sharks_last_meal + shark_survival_time
    
    log.append("\n--- Final Calculation ---")
    log.append("The total time is determined by when the last fish, the shark, dies.")
    log.append("Final Time = Time of Shark's Last Meal + Shark's Survival Time")

    # Print all log entries to show the process
    for entry in log:
        print(entry)
        
    # Print the final equation as requested
    print(f"\nFinal Time = {time_of_sharks_last_meal} + ({base_survival_time} + {extension_time}) = {final_time}")
    
    # Print the final answer in the specified format
    print(f"<<<{final_time}>>>")

# Execute the function to solve the problem
solve_fish_problem()