def solve_fish_problem():
    """
    This function simulates the optimal strategy to clear the pond of all fish
    in the shortest possible time and prints the step-by-step logic.
    """
    # Initial State
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eaten_count = 0
    time = 0

    print(f"Time: {time} minutes")
    print(f"Initial state: {goldfish_in_pond} goldfish in the pond, 1 shark. The fisherman's basket is empty.")
    print("-" * 30)

    # --- Phase 1: t=0 to t=15 ---
    # The fisherman uses the fast Rod B to catch 3 fish.
    # The shark eats on its own schedule (every 10 minutes).
    print("Phase 1: Fisherman uses the fast Rod B.")
    
    # t=5: Fisherman catches 1st fish
    time_event_1 = 5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    
    # t=10: Fisherman catches 2nd fish AND shark eats 1 from pond
    time_event_2 = 10
    goldfish_in_pond -= 1  # Fisherman
    goldfish_in_basket += 1
    goldfish_in_pond -= 1  # Shark
    shark_eaten_count += 1
    
    # t=15: Fisherman catches 3rd fish
    time = 15
    goldfish_in_pond -= 1
    goldfish_in_basket += 1

    print(f"Time: {time} minutes")
    print("The fisherman has caught 3 goldfish. The shark has eaten 1.")
    print(f"State: {goldfish_in_pond} goldfish in pond, {goldfish_in_basket} in basket, shark has eaten {shark_eaten_count} total.")
    print("-" * 30)

    # --- Phase 2: t=15 to t=22 ---
    # To avoid the basket overflowing, the fisherman switches to Rod A.
    # This gives him a free hand to feed the 3 caught fish to the shark.
    # This takes 3 * 2 = 6 minutes of feeding time (t=15 to t=21).
    # Using Rod A for 7 minutes (t=15 to t=22) also resets Rod B's usage timer.
    print("Phase 2: Fisherman switches to Rod A to feed the shark and manage the basket.")
    
    goldfish_in_basket -= 3
    shark_eaten_count += 3
    
    # t=22: Fisherman catches 1 fish with Rod A.
    time = 22
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    
    print(f"Time: {time} minutes")
    print("The fisherman fed 3 fish to the shark and caught 1 more with Rod A.")
    print(f"State: {goldfish_in_pond} goldfish in pond, {goldfish_in_basket} in basket, shark has eaten {shark_eaten_count} total.")
    print("-" * 30)

    # --- Phase 3: t=22 to t=32 ---
    # Fisherman switches back to the faster Rod B.
    print("Phase 3: Fisherman switches back to Rod B.")

    # t=27: Fisherman catches a fish
    time_event_3 = 27
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    
    # t=31: Shark eats 1 from pond (10 mins after its last meal at t=21, from feeding)
    time_event_4 = 31
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    
    # t=32: Fisherman catches another fish
    time = 32
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    
    print(f"Time: {time} minutes")
    print("The fisherman caught 2 more fish with Rod B, and the shark ate 1 from the pond.")
    print(f"State: {goldfish_in_pond} goldfish in pond, {goldfish_in_basket} in basket, shark has eaten {shark_eaten_count} total.")
    print("-" * 30)

    # --- Phase 4: t=32 to t=38 ---
    # The endgame. The key is to upgrade the shark.
    print("Phase 4: Upgrading the shark and clearing the pond.")
    
    # t=32-34: Fisherman switches to Rod A and feeds one fish.
    time_event_5 = 34
    goldfish_in_basket -= 1
    shark_eaten_count += 1
    print(f"Time: {time_event_5} minutes")
    print(f"Fisherman feeds the shark its 6th meal. The shark is now stronger and will eat every 2 minutes.")
    
    # Now, let the strong shark clear the remaining 2 fish from the pond.
    # t=36: Shark eats 1 from pond
    time_event_6 = 36
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    
    # t=38: Shark eats the final fish from the pond.
    time_pond_empty = 38
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    
    print(f"Time: {time_pond_empty} minutes")
    print("The upgraded shark eats the last two goldfish from the pond.")
    print(f"State: All {goldfish_in_pond} goldfish in the pond are gone. The shark has eaten {shark_eaten_count} total.")
    print("-" * 30)

    # --- Phase 5: Shark Starvation ---
    # The shark's last meal was at t=38. It has eaten 8 fish in the last 60 minutes (>4).
    # Its survival time is extended by 4 minutes, from 11 to 15.
    print("Phase 5: Waiting for the shark to starve.")
    base_starvation_time = 11
    extended_starvation_bonus = 4
    starvation_time = base_starvation_time + extended_starvation_bonus
    
    print(f"The shark's last meal was at {time_pond_empty} minutes.")
    print(f"Because it ate more than 4 fish in the last hour, its survival time is {starvation_time} minutes.")
    
    final_time = time_pond_empty + starvation_time
    
    print("\nFinal Calculation:")
    print(f"The earliest time for no fish left is the time the pond is empty plus the shark's starvation time.")
    print(f"Final Time = Time Pond Empty + Shark Starvation Time")
    # The final equation with each number printed out
    print(f"{time_pond_empty} + {starvation_time} = {final_time}")
    
solve_fish_problem()