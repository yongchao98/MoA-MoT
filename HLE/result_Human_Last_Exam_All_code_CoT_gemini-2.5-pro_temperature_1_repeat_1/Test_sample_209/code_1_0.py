def solve_fish_pond_problem():
    """
    Calculates the earliest time the pond is empty by simulating the optimal strategy.
    """
    # Initial State
    initial_goldfish = 10
    time = 0
    
    # State variables
    goldfish_in_pond = initial_goldfish
    fisherman_caught_total = 0
    shark_eaten_total = 0
    
    print("Starting simulation...")
    print(f"Time {time}: The pond has {goldfish_in_pond} goldfish and 1 shark.")
    print("-" * 30)

    # --- Phase 1: Fisherman uses Rod B until the basket is full ---
    print("Phase 1: Fisherman uses Rod B (5 min/fish) to fill his basket (capacity 4).")
    
    catches_to_fill_basket = 4
    time_to_fill_basket = catches_to_fill_basket * 5  # 4 fish * 5 min/fish
    
    # During this time, the shark eats every 10 minutes.
    shark_eats_p1 = time_to_fill_basket // 10
    
    # Update state at the end of Phase 1
    time += time_to_fill_basket
    fisherman_caught_p1 = catches_to_fill_basket
    fisherman_caught_total += fisherman_caught_p1
    shark_eaten_total += shark_eats_p1
    
    # Final equation for Phase 1
    goldfish_removed_p1 = fisherman_caught_p1 + shark_eats_p1
    goldfish_in_pond -= goldfish_removed_p1
    
    print(f"Fisherman catches {fisherman_caught_p1} fish in {time_to_fill_basket} minutes (at t=5, 10, 15, 20).")
    print(f"During this time, the shark eats {shark_eats_p1} fish (at t=10, 20).")
    print(f"Time {time}: Total fish removed so far = {fisherman_caught_p1} (caught) + {shark_eats_p1} (eaten) = {goldfish_removed_p1}.")
    print(f"           Pond now has {initial_goldfish} - {goldfish_removed_p1} = {goldfish_in_pond} goldfish.")
    print(f"           The basket is full. Fisherman leaves to get a new one.")
    print("-" * 30)

    # --- Phase 2: Fisherman is away getting a new basket ---
    print("Phase 2: Fisherman is away for 20 minutes to get a new basket.")
    
    time_away = 20
    
    # The shark's last meal was at t=20. It eats every 10 minutes.
    shark_eats_p2 = time_away // 10
    
    # Update state at the end of Phase 2
    time += time_away
    shark_eaten_total += shark_eats_p2
    
    # Final equation for Phase 2
    goldfish_removed_p2 = shark_eats_p2
    goldfish_in_pond -= goldfish_removed_p2
    
    print(f"During the {time_away} minutes the fisherman is away, the shark eats {shark_eats_p2} more fish (at t=30, 40).")
    print(f"Time {time}: Fisherman returns. The pond now has {goldfish_in_pond + goldfish_removed_p2} - {goldfish_removed_p2} = {goldfish_in_pond} goldfish.")
    print("-" * 30)

    # --- Phase 3: Fisherman fishes for the remaining goldfish ---
    print(f"Phase 3: Fisherman uses Rod B to catch the remaining {goldfish_in_pond} fish.")
    
    remaining_fish = goldfish_in_pond
    time_to_catch_remaining = remaining_fish * 5 # Using Rod B
    
    # During this final fishing session (from t=40 to t=50), the shark's next meal is due at t=50.
    # At t=50, the fisherman makes his final catch, emptying the pond.
    # We assume the fisherman, acting optimally, catches the fish at the same moment the shark would have.
    shark_eats_p3 = 0
    fisherman_caught_p3 = remaining_fish
    
    # Update state at the end of Phase 3
    time += time_to_catch_remaining
    fisherman_caught_total += fisherman_caught_p3
    
    # Final equation for Phase 3
    goldfish_removed_p3 = fisherman_caught_p3
    goldfish_in_pond -= goldfish_removed_p3

    print(f"It takes the fisherman {time_to_catch_remaining} more minutes to catch the last {remaining_fish} fish (at t=45, 50).")
    print(f"Time {time}: The last goldfish is caught. The pond is empty.")
    print("-" * 30)
    
    # --- Final Summary ---
    print("Final Tally:")
    print(f"Total fish caught by fisherman: {fisherman_caught_total}")
    print(f"Total fish eaten by shark: {shark_eaten_total}")
    print(f"Total goldfish removed from pond: {fisherman_caught_total} + {shark_eaten_total} = {initial_goldfish}")
    
    print("\nThe earliest possible time when there are no more fish left in the pond is:")
    print(f"{time} minutes")

solve_fish_pond_problem()