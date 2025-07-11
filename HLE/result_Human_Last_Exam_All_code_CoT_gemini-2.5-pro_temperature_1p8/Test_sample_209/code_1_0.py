def solve_fish_pond_problem():
    """
    Simulates the optimal strategy to clear the pond of goldfish and prints the timeline.
    """
    initial_goldfish = 10
    
    print("### Finding the Earliest Time to Clear the Pond ###\n")
    print("The optimal strategy is for the fisherman to use the faster Rod B to catch fish as quickly as possible.")
    print("This will fill the basket, forcing a 20-minute trip for a new one, but this is faster than using the slower rod.\n")
    
    print("--- Timeline of Events ---")
    
    # Phase 1: Fishing until the basket is full
    print("\n[Phase 1: Fisherman uses Rod B until basket is full (t=0 to t=20)]")
    # Fisherman catches at t=5, 10, 15, 20
    fisherman_catches_phase1 = 4
    print(f"Time 5, 10, 15, 20: The fisherman catches {fisherman_catches_phase1} goldfish.")
    # Shark eats at t=10, 20
    shark_eats_phase1 = 2
    print(f"Time 10, 20: The shark eats {shark_eats_phase1} goldfish from the pond.")
    
    goldfish_at_t20 = initial_goldfish - fisherman_catches_phase1 - shark_eats_phase1
    print(f"By Time 20, the fisherman's basket is full (4 fish).")
    print(f"Number of goldfish left in the pond: {initial_goldfish} - {fisherman_catches_phase1} - {shark_eats_phase1} = {goldfish_at_t20}")
    
    # Phase 2: Fisherman's trip
    print("\n[Phase 2: Fisherman's 20-minute basket trip (t=20 to t=40)]")
    print("The fisherman is away. Meanwhile, the shark continues to eat every 10 minutes.")
    # Shark eats at t=30, 40
    shark_eats_phase2 = 2
    print("Time 30: The shark eats 1 goldfish.")
    print("Time 40: The shark eats 1 goldfish.")
    
    goldfish_at_t40 = goldfish_at_t20 - shark_eats_phase2
    print(f"When the fisherman returns at Time 40, the number of goldfish is: {goldfish_at_t20} - {shark_eats_phase2} = {goldfish_at_t40}")

    # Phase 3: Final catches
    print("\n[Phase 3: Clearing the remaining fish (t=40 onwards)]")
    print("The fisherman returns with an empty basket and resumes fishing with Rod B.")
    fisherman_catches_phase3 = 1
    print("Time 45: The fisherman catches 1 goldfish.")
    goldfish_after_catch = goldfish_at_t40 - fisherman_catches_phase3
    print(f"Goldfish remaining: {goldfish_at_t40} - {fisherman_catches_phase3} = {goldfish_after_catch}")

    shark_eats_phase3 = 1
    final_time = 50
    print(f"Time {final_time}: The shark's 10-minute cycle is up, and it eats the last remaining goldfish.")
    
    final_goldfish = goldfish_after_catch - shark_eats_phase3
    
    print("\n--- Final Calculation ---")
    total_fisherman_catches = fisherman_catches_phase1 + fisherman_catches_phase3
    total_shark_eats = shark_eats_phase1 + shark_eats_phase2 + shark_eats_phase3
    
    print("The final equation for the number of goldfish removed from the pond is:")
    print(f"Initial Goldfish - Goldfish Caught by Fisherman - Goldfish Eaten by Shark = Final Goldfish")
    print(f"{initial_goldfish} - {total_fisherman_catches} - {total_shark_eats} = {final_goldfish}")
    
    print(f"\nThe last goldfish was removed from the pond at Time {final_time}.")

solve_fish_pond_problem()
print("\n<<<50>>>")