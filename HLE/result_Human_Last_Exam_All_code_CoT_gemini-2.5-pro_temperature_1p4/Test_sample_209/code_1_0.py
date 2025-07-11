def solve_fish_problem():
    """
    Simulates the optimal strategy to find the earliest time the pond is empty.
    The optimal strategy is to use the fastest fishing rod (Rod B) and incur
    the 20-minute basket-replacement penalty, as this is faster than
    slowing down to feed the shark.
    """
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eats_count = 0
    shark_last_meal_time = 0

    print(f"Time {time}: The simulation begins. There are {goldfish_in_pond} goldfish and 1 shark in the pond.")
    print("-" * 20)

    # Fisherman uses Rod B (5 min/fish) and shark eats on its own (10 min/fish)
    print("Phase 1: Fisherman uses Rod B. Shark eats on its own schedule.")

    # Time 5
    time = 5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish with Rod B. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")

    # Time 10 (Simultaneous events)
    time = 10
    goldfish_in_pond -= 1  # Shark eats
    shark_eats_count += 1
    shark_last_meal_time = 10
    print(f"Time {time}: Shark eats a free-swimming goldfish. Pond now: {goldfish_in_pond}. Shark has eaten: {shark_eats_count}.")
    goldfish_in_pond -= 1  # Fisherman catches
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish with Rod B. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")

    # Time 15
    time = 15
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish with Rod B. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")

    # Time 20 (Simultaneous events)
    time = 20
    goldfish_in_pond -= 1  # Shark eats
    shark_eats_count += 1
    shark_last_meal_time = 20
    print(f"Time {time}: Shark eats a free-swimming goldfish. Pond now: {goldfish_in_pond}. Shark has eaten: {shark_eats_count}.")
    goldfish_in_pond -= 1  # Fisherman catches
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches a goldfish with Rod B. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}. The basket is now full.")
    print("-" * 20)
    
    # Phase 2: Basket penalty trip
    print("Phase 2: Basket overflow triggers a 20-minute trip for a new basket.")

    # Time 25
    time = 25
    goldfish_in_pond -= 1
    print(f"Time {time}: Fisherman catches a goldfish, overflowing the basket. Pond: {goldfish_in_pond}.")
    print(f"Time {time}: He starts the 20-minute trip and will return at Time {time + 20}.")
    goldfish_in_basket = 0  # He gets a new, empty basket.

    # Time 30
    time = 30
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_meal_time = 30
    print(f"Time {time}: While fisherman is away, the shark eats a goldfish. Pond: {goldfish_in_pond}. Shark has eaten: {shark_eats_count}.")

    # Time 40
    time = 40
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_meal_time = 40
    print(f"Time {time}: The shark eats another goldfish. Pond: {goldfish_in_pond}. Shark has eaten: {shark_eats_count}.")
    print("-" * 20)

    # Phase 3: Final catch and shark's demise
    print("Phase 3: Fisherman returns to catch the final goldfish.")
    
    # Time 45
    time = 45
    print(f"Time {time}: The fisherman returns with a new basket.")

    # Time 50
    time = 50
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time}: Fisherman catches the last goldfish. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")
    print("All goldfish are now gone from the pond.")
    print("-" * 20)
    
    # Final Calculation
    print("Final Calculation: Determining the shark's death.")
    starvation_time = 11
    # Check for survival extension: eaten > 4 fish.
    # The shark ate exactly 4 fish, so the condition is not met.
    print(f"The shark's last meal was at Time {shark_last_meal_time}.")
    print(f"It ate a total of {shark_eats_count} fish. This does not trigger the survival extension (> 4 needed).")
    print(f"The standard starvation time is {starvation_time} minutes.")

    shark_death_time = shark_last_meal_time + starvation_time
    print(f"The shark will die at: {shark_last_meal_time} (last meal) + {starvation_time} (starvation) = {shark_death_time} minutes.")
    
    print("\n---------------------------------------------------------------------------------")
    print(f"The earliest possible time when there are no more fish left in the pond is {shark_death_time} minutes.")
    print("---------------------------------------------------------------------------------")

solve_fish_problem()