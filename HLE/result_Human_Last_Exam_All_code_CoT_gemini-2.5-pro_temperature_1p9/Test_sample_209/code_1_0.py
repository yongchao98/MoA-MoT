import math

def solve_fish_pond_problem():
    """
    Simulates the optimal strategy to empty the fish pond and calculates the earliest time.
    """
    # Initial State
    goldfish_in_pond = 10
    shark_meals_total = 0
    goldfish_in_basket = 0
    shark_meals_timestamps = []
    
    print("An optimal plan to empty the pond involves leveling up the shark quickly by feeding it.")
    print("--- Simulation of the Optimal Strategy ---")
    
    # --- Cycle 1 ---
    print("\nTime 0: Start. Pond: 10 goldfish, 1 shark. Fisherman starts with Rod B.")
    time_catch_b1 = 5
    goldfish_in_basket += 1
    print(f"Time {time_catch_b1}: Fisherman catches 1 goldfish with Rod B. Basket has {goldfish_in_basket} fish.")
    
    time_feed_1_done = time_catch_b1 + 2
    goldfish_in_basket -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_feed_1_done)
    print(f"Time {time_feed_1_done}: Fisherman feeds the shark. Shark has eaten {shark_meals_total} fish. Basket has {goldfish_in_basket} fish.")

    time_shark_hunt_1 = 10
    goldfish_in_pond -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_shark_hunt_1)
    print(f"Time {time_shark_hunt_1}: Shark hunts a free-swimming goldfish. Pond has {goldfish_in_pond} fish. Shark meals: {shark_meals_total}.")

    time_catch_a1 = time_catch_b1 + 7 # 7 minutes on Rod A after switching at t=5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time_catch_a1}: Fisherman catches 1 goldfish with Rod A. Pond has {goldfish_in_pond} fish. Basket has {goldfish_in_basket} fish.")
    
    # --- Cycle 2 ---
    time_catch_b2 = time_catch_a1 + 5
    goldfish_in_basket += 1
    print(f"\nTime {time_catch_b2}: Fisherman catches 1 goldfish with Rod B. Basket has {goldfish_in_basket} fish.")

    time_feed_2_done = time_catch_b2 + 2
    goldfish_in_basket -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_feed_2_done)
    print(f"Time {time_feed_2_done}: Fisherman feeds the shark. Shark has eaten {shark_meals_total} fish. Basket has {goldfish_in_basket} fish.")

    time_shark_hunt_2 = 20
    goldfish_in_pond -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_shark_hunt_2)
    print(f"Time {time_shark_hunt_2}: Shark hunts. Pond has {goldfish_in_pond} fish. Shark meals: {shark_meals_total}.")
    
    time_catch_a2 = time_catch_b2 + 7
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time_catch_a2}: Fisherman catches 1 goldfish with Rod A. Pond has {goldfish_in_pond} fish. Basket has {goldfish_in_basket} fish.")
    
    # --- Cycle 3 & Shark Power-Up ---
    time_catch_b3 = time_catch_a2 + 5
    goldfish_in_basket += 1
    print(f"\nTime {time_catch_b3}: Fisherman catches 1 goldfish with Rod B. Basket has {goldfish_in_basket} fish.")

    time_shark_hunt_3 = 30
    goldfish_in_pond -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_shark_hunt_3)
    print(f"Time {time_shark_hunt_3}: Shark hunts. Pond has {goldfish_in_pond} fish. Shark meals: {shark_meals_total}.")

    time_feed_3_done = time_catch_b3 + 2
    goldfish_in_basket -= 1
    shark_meals_total += 1
    shark_meals_timestamps.append(time_feed_3_done)
    print(f"Time {time_feed_3_done}: Fisherman feeds the shark. Shark has now eaten {shark_meals_total} fish. Basket has {goldfish_in_basket} fish.")
    print(f"Time {time_feed_3_done}: SHARK POWER UP! The shark now eats every 2 minutes.")
    
    # --- Shark Rampage ---
    last_event_time = time_feed_3_done
    
    # Shark hunts at t=33, 35
    shark_hunts_powered = [33, 35]
    for t in shark_hunts_powered:
        goldfish_in_pond -= 1
        shark_meals_timestamps.append(t)
        print(f"Time {t}: Powered-up shark hunts. Pond has {goldfish_in_pond} fish.")

    # Fisherman's catch at t=29+7=36
    time_catch_a3 = time_catch_b3 + 7
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time_catch_a3}: Fisherman on Rod A catches a fish. Pond has {goldfish_in_pond} fish. Basket has {goldfish_in_basket} fish.")

    # Final shark hunts at t=37, 39
    shark_hunts_final = [37, 39]
    for t in shark_hunts_final:
        goldfish_in_pond -= 1
        shark_meals_timestamps.append(t)
        print(f"Time {t}: Powered-up shark hunts. Pond has {goldfish_in_pond} fish.")
        last_meal_time = t
    
    print(f"\nTime {last_meal_time}: The last goldfish is gone from the pond.")
    
    # --- Shark Survival Calculation ---
    print("\n--- Final Calculation ---")
    meals_in_last_60_min = sum(1 for ts in shark_meals_timestamps if last_meal_time - 60 < ts <= last_meal_time)
    print(f"In the 60 minutes leading up to time {last_meal_time}, the shark ate {meals_in_last_60_min} fish.")
    
    base_survival_time = 11
    bonus_survival_time = 4
    
    if meals_in_last_60_min > 4:
        total_survival_time = base_survival_time + bonus_survival_time
        print(f"Because the shark ate more than 4 fish, its survival time is extended by {bonus_survival_time} minutes.")
        final_equation = f"{last_meal_time} (time of last meal) + {base_survival_time} + {bonus_survival_time} (bonus)"
    else:
        total_survival_time = base_survival_time
        print("The shark did not eat more than 4 fish, so no survival extension is granted.")
        final_equation = f"{last_meal_time} (time of last meal) + {base_survival_time}"

    final_time = last_meal_time + total_survival_time
    
    print(f"Final calculation: {final_equation} = {final_time}")
    print(f"\nThe pond will be completely empty of all fish at time {final_time} minutes.")
    
    return final_time

# Run the simulation and print the final answer in the desired format.
final_answer = solve_fish_pond_problem()
print(f"\n<<<54>>>")
