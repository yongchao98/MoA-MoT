def solve_fish_pond_problem():
    """
    Calculates the earliest time all fish are gone from the pond by simulating an optimal strategy.
    """
    # Initial state
    time = 0
    pond_gf = 10
    basket_gf = 0
    shark_eaten_count = 0
    # Store tuples of (time, source) for shark's meals
    shark_eaten_log = [] 
    shark_starvation_period = 11
    
    print("--- Strategy: Power-up the shark by feeding it while managing resources efficiently ---\n")

    # === Phase 1: Time 0-15 ===
    # Fisherman uses Rod B to quickly catch 3 fish.
    print(f"Time 0: Start. Pond has {pond_gf} goldfish. Fisherman uses Rod B.")
    
    # t=5
    time = 5
    pond_gf -= 1
    basket_gf += 1
    print(f"Time {time}: Fisherman (Rod B) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")
    
    # t=10
    time = 10
    # Fisherman's catch
    pond_gf -= 1 
    basket_gf += 1
    print(f"Time {time}: Fisherman (Rod B) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")
    # Shark's natural meal
    pond_gf -= 1 
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Shark eats a free-swimming goldfish. Pond={pond_gf}, Shark has eaten {shark_eaten_count}.")

    # t=15
    time = 15
    pond_gf -= 1
    basket_gf += 1
    print(f"Time {time}: Fisherman (Rod B) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")
    print(f"INFO: Basket has 3 fish. Fisherman switches to Rod A to feed the shark.\n")

    # === Phase 2: Time 15-22 ===
    # Fisherman uses Rod A. He can now feed the shark from his basket.
    # Feeding (2 mins) can overlap with fishing with Rod A (7 mins).
    
    # t=17
    time = 17
    basket_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Fisherman feeds a caught fish to the shark. Basket={basket_gf}, Shark has eaten {shark_eaten_count}.")
    
    # t=19
    time = 19
    basket_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Fisherman feeds a caught fish to the shark. Basket={basket_gf}, Shark has eaten {shark_eaten_count}.")
    
    # t=21
    time = 21
    basket_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Fisherman feeds a caught fish to the shark. Basket={basket_gf}, Shark has eaten {shark_eaten_count}.")

    # t=22
    time = 22 
    pond_gf -= 1
    basket_gf += 1
    print(f"Time {time}: After 7 minutes, Fisherman (Rod A) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")
    print(f"INFO: Fisherman switches back to the faster Rod B.\n")

    # === Phase 3: Time 22-38 ===
    # Fisherman on Rod B, shark getting closer to power-up.
    
    # t=27
    time = 27
    pond_gf -= 1
    basket_gf += 1
    print(f"Time {time}: Fisherman (Rod B) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")

    # t=31
    time = 31 
    pond_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Shark eats a free-swimming goldfish. Pond={pond_gf}, Shark has eaten {shark_eaten_count}.")
    
    # Check for starvation extension rule
    eaten_in_last_60 = len([t for t in shark_eaten_log if time - t <= 60])
    if eaten_in_last_60 > 4 and shark_starvation_period == 11:
        shark_starvation_period = 15
        print(f"INFO: Shark has now eaten {eaten_in_last_60} fish. Starvation time is extended to {shark_starvation_period} minutes.")

    # t=32
    time = 32
    pond_gf -= 1
    basket_gf += 1
    print(f"Time {time}: Fisherman (Rod B) catches a goldfish. Pond={pond_gf}, Basket={basket_gf}.")
    print(f"INFO: Basket has 3 fish. Fisherman switches to Rod A.\n")

    # t=34: Shark Power-up!
    time = 34
    basket_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    print(f"Time {time}: Fisherman feeds a caught fish to the shark. Basket={basket_gf}, Shark has eaten {shark_eaten_count}.")
    print(f"**POWER UP**: Shark has eaten its 6th fish! It now eats every 2 minutes.")

    # t=36: Powered-up shark eats
    time = 36
    pond_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    last_shark_meal_time = time
    print(f"Time {time}: Powered-up shark eats a free-swimming goldfish. Pond={pond_gf}.")

    # t=38: Powered-up shark eats the last goldfish
    time = 38
    pond_gf -= 1
    shark_eaten_count += 1
    shark_eaten_log.append(time)
    last_shark_meal_time = time
    print(f"Time {time}: Powered-up shark eats the last free-swimming goldfish. Pond={pond_gf}.")
    
    print("\n--- All 10 goldfish are gone from the pond. ---")
    
    # === Final Phase: Shark Starvation ===
    print("The final step is to wait for the shark to starve.")
    
    final_time = last_shark_meal_time + shark_starvation_period
    
    print(f"\nThe shark's last meal was at time {last_shark_meal_time} minutes.")
    print(f"Because it ate >4 fish in the last hour, its survival time is {shark_starvation_period} minutes.")
    print("\nThe final calculation is:")
    print(f"Time of last meal + Shark starvation period = Total time")
    print(f"{last_shark_meal_time} minutes + {shark_starvation_period} minutes = {final_time} minutes")
    
    return final_time

# Run the simulation and print the final answer
final_answer = solve_fish_pond_problem()
print(f"\n<<<The earliest possible time when there are no more fish left in the pond is {final_answer} minutes.>>>")
