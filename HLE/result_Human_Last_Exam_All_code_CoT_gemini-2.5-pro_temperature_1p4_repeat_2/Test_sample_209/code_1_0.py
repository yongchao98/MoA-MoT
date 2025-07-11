def solve_fish_problem():
    """
    This function simulates the optimal strategy for the fisherman to empty the pond
    and prints a step-by-step timeline of the events.
    """
    # Initial State
    time = 0
    free_fish = 10
    basket_fish = 0
    shark_hunted_fish = 0
    shark_fed_fish = 0

    def print_status(event_time, description):
        total_shark_eaten = shark_hunted_fish + shark_fed_fish
        print(f"--- Time: {event_time} minutes ---")
        print(description)
        print(f"State: {free_fish} free fish, {basket_fish} in basket, {total_shark_eaten} eaten by shark ({shark_hunted_fish} hunted, {shark_fed_fish} fed).\n")

    print_status(time, "The fisherman begins. He chooses the fastest rod, Rod B, to catch fish as quickly as possible.")

    # Step 1: Fisherman catches first fish with Rod B
    time += 5  # Time for Rod B
    free_fish -= 1
    basket_fish += 1
    print_status(time, "The fisherman catches his first fish with Rod B (5 minutes).")

    # Step 2: Fisherman catches second fish, shark eats its first
    time += 5  # 10 minutes total
    free_fish -= 1  # Fisherman catches
    basket_fish += 1
    free_fish -= 1  # Shark hunts (happens at the 10-minute mark)
    shark_hunted_fish += 1
    print_status(time, "The fisherman catches his second fish (at t=10). Simultaneously, the shark hunts its first free-swimming fish (every 10 minutes).")

    # Step 3: Fisherman catches third fish
    time += 5  # 15 minutes total
    free_fish -= 1
    basket_fish += 1
    print_status(time, "The fisherman catches his third fish (at t=15). The basket now has 3 fish. The next catch would make it 4, forcing a 20-minute trip. This must be avoided.")

    # Step 4: Switch to Rod A to feed the shark and manage basket
    time += 2  # 17 minutes total
    basket_fish -= 1
    shark_fed_fish += 1
    description = (
        "Strategy: To free up basket space, the fisherman switches to Rod A (one-handed) and feeds a caught fish to the shark.\n"
        "This takes 2 minutes. This action resets the shark's 10-minute hunting timer."
    )
    print_status(time, description)

    # Step 5: Wait to reset Rod B timer and switch back
    time += 5 # 22 minutes total (15 + 2 + 5 = 22, total of 7 mins on Rod A)
    description = (
        "To use the fast Rod B again, the fisherman must hold Rod A for 7 minutes. He started at t=15, so he can switch back at t=22.\n"
        "He now switches back to Rod B."
    )
    print_status(time, description)

    # Step 6: Catch another fish and shark hunts again
    time += 5 # 27 minutes total
    free_fish -= 1
    basket_fish += 1
    free_fish -= 1 # Shark's next hunt is at t=17 + 10 = 27
    shark_hunted_fish += 1
    print_status(time, "Fisherman catches a fish with Rod B. Simultaneously, the shark hunts again. The basket is back at 3 fish.")

    # Step 7: Repeat the feed cycle
    time += 2 # 29 minutes total
    basket_fish -= 1
    shark_fed_fish += 1
    print_status(time, "Again, the fisherman switches to Rod A and feeds the shark to manage basket space (2 mins). This resets the shark's hunt timer again.")

    # Step 8: Wait and switch back to Rod B
    time += 5 # 34 minutes total (27 + 2 + 5 = 34)
    print_status(time, "The fisherman holds Rod A for the required 7 minutes (from t=27 to t=34) and switches back to Rod B.")

    # Step 9: Catch another fish and shark hunts
    time += 5 # 39 minutes total
    free_fish -= 1
    basket_fish += 1
    free_fish -= 1 # Shark's next hunt is at t=29 + 10 = 39
    shark_hunted_fish += 1
    total_shark_eaten = shark_hunted_fish + shark_fed_fish
    description = (
        f"The fisherman catches another fish (now has {basket_fish} in basket). The shark also hunts, having now eaten a total of {total_shark_eaten} fish.\n"
        f"There are now only {free_fish} fish left in the pond."
    )
    print_status(time, description)
    
    # Step 10: Final crucial decision - Power up the shark
    time += 2 # 41 minutes total
    basket_fish -= 1
    shark_fed_fish += 1
    total_shark_eaten = shark_hunted_fish + shark_fed_fish
    description = (
        "CRUCIAL MOVE: The shark has eaten 5 fish. The fisherman can feed it a 6th fish to activate its power-up.\n"
        f"The fisherman feeds the shark (2 mins). The shark has now eaten {total_shark_eaten} fish and will hunt every 2 minutes."
    )
    print_status(time, description)
    
    # Step 11: Super-shark cleans up
    time += 2 # 43 minutes total
    free_fish -= 1
    shark_hunted_fish += 1
    print_status(time, f"The powered-up shark eats a fish (after 2 minutes). Only {free_fish} fish remains.")
    
    # Step 12: Shark eats the final fish
    time += 2 # 45 minutes total
    free_fish -= 1
    shark_hunted_fish += 1
    print_status(time, f"The shark eats the last fish in the pond (another 2 minutes). The pond is now empty.")
    
    print("--------------------------------------------------")
    print(f"The earliest possible time when there are no more fish left in the pond is {time} minutes.")
    print("--------------------------------------------------")

solve_fish_problem()