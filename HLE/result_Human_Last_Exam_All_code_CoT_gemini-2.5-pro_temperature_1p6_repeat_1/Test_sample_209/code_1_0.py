def solve_fish_pond_problem():
    """
    Simulates the optimal strategy to find the earliest time the pond is empty.
    """
    log = []
    
    # State variables
    time = 0
    goldfish_pond = 10
    basket = 0
    shark_eaten = 0
    shark_last_meal_time = 0

    log.append(f"T={time} min: Start. The pond has 10 goldfish and 1 shark.")
    log.append("Fisherman's optimal strategy is to use the fastest tool, Rod B (5 min/fish).")
    
    # --- Phase 1: Fishing until the basket is full ---
    log.append("\n--- Phase 1: Fishing with Rod B ---")
    
    # Event at T=5
    time = 5
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches a fish. Pond: {goldfish_pond}, Basket: {basket}.")
    
    # Event at T=10 (simultaneous)
    time = 10
    # Fisherman's catch
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches a fish. Pond is now {goldfish_pond}, Basket is {basket}.")
    # Shark's meal
    goldfish_pond -= 1
    shark_eaten += 1
    shark_last_meal_time = time
    log.append(f"T={time} min: Shark eats a fish. Pond is now {goldfish_pond}. Shark has eaten {shark_eaten} fish.")

    # Event at T=15
    time = 15
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches a fish. Pond: {goldfish_pond}, Basket: {basket}.")
    
    # Event at T=20 (simultaneous)
    time = 20
    # Fisherman's catch
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches a fish. Pond is now {goldfish_pond}, Basket is {basket}.")
    # Shark's meal
    goldfish_pond -= 1
    shark_eaten += 1
    shark_last_meal_time = time
    log.append(f"T={time} min: Shark eats a fish. Pond is now {goldfish_pond}. Shark has eaten {shark_eaten} fish.")
    
    # --- Basket Trip ---
    log.append("\n--- Phase 2: Basket Trip ---")
    log.append(f"T={time} min: The basket is full ({basket}/4). Fisherman starts a 20-minute trip to get a new one.")
    trip_start_time = time
    
    # Shark's actions during the trip
    time_of_next_shark_meal = shark_last_meal_time + 10
    time = time_of_next_shark_meal
    goldfish_pond -= 1
    shark_eaten += 1
    shark_last_meal_time = time
    log.append(f"T={time} min: While fisherman is away, shark eats a fish. Pond: {goldfish_pond}.")

    time_of_next_shark_meal = shark_last_meal_time + 10
    time = time_of_next_shark_meal
    goldfish_pond -= 1
    shark_eaten += 1
    shark_last_meal_time = time
    log.append(f"T={time} min: While fisherman is away, shark eats another fish. Pond: {goldfish_pond}.")
    
    # Fisherman returns
    time = trip_start_time + 20
    basket = 0
    log.append(f"T={time} min: Fisherman returns with an empty basket.")

    # --- Phase 3: Catching the remaining fish ---
    log.append("\n--- Phase 3: Final Catches ---")
    
    # Next catch
    time = time + 5
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches a fish. Pond: {goldfish_pond}, Basket: {basket}.")

    # Final catch
    time = time + 5
    goldfish_pond -= 1
    basket += 1
    log.append(f"T={time} min: Fisherman catches the last goldfish. The pond is now empty of goldfish.")
    
    # --- Final Calculation ---
    log.append("\n--- Final Calculation: Emptying the Pond ---")
    # Check shark's starvation rule. It starves if it doesn't eat for 11 minutes.
    # The rule for extended survival (more than 4 fish in 60 mins) is not met, as it ate exactly 4.
    starvation_duration = 11
    
    log.append(f"The last time the shark ate was at T={shark_last_meal_time} minutes.")
    log.append(f"With no goldfish left, the shark cannot eat and will starve in {starvation_duration} minutes.")
    
    final_time = shark_last_meal_time + starvation_duration
    
    # Print the full story
    for event in log:
        print(event)
        
    print("\n=======================================================")
    print("The final calculation to determine when no fish are left:")
    print(f"Time of shark's last meal: {shark_last_meal_time} minutes")
    print(f"Time until shark starves: {starvation_duration} minutes")
    print(f"Final Time = Last Meal Time + Starvation Duration")
    print(f"Final Time = {shark_last_meal_time} + {starvation_duration} = {final_time}")
    print("=======================================================")

if __name__ == '__main__':
    solve_fish_pond_problem()
    print("\n<<<51>>>")