def solve_fish_pond_problem():
    """
    Calculates and explains the earliest possible time when there are no more fish left in the pond.
    """
    # Initial state
    goldfish_in_pond = 10
    shark_eaten_count = 0
    basket_count = 0
    
    print("Finding the earliest time for the pond to be empty.")
    print("The optimal strategy is for the fisherman to use the fastest rod (Rod B) as much as possible.")
    print("---")
    
    # Step-by-step simulation of the optimal strategy
    print(f"Time 0: The simulation starts. There are {goldfish_in_pond} goldfish in the pond.")
    print("The fisherman starts using Rod B (catches a fish every 5 minutes).")
    print("The shark eats a free-swimming goldfish every 10 minutes (at t=10, 20, 30...).")
    print("---")

    # Time 5
    time = 5
    goldfish_in_pond -= 1
    basket_count += 1
    print(f"Time {time}: The fisherman catches a goldfish. {goldfish_in_pond} goldfish are left in the pond, and {basket_count} is in the basket.")
    print("---")

    # Time 10
    time = 10
    # Fisherman's catch
    goldfish_in_pond -= 1
    basket_count += 1
    # Shark's meal
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 10
    print(f"Time {time}: The fisherman catches another goldfish, and the shark eats one.")
    print(f"There are now {goldfish_in_pond} goldfish in the pond, and {basket_count} in the basket.")
    print("---")

    # Time 15
    time = 15
    goldfish_in_pond -= 1
    basket_count += 1
    print(f"Time {time}: The fisherman catches a goldfish. {goldfish_in_pond} goldfish are left in the pond, and {basket_count} are in the basket.")
    print("---")

    # Time 20
    time = 20
    # Fisherman's catch
    goldfish_in_pond -= 1
    basket_count += 1
    # Shark's meal
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 20
    print(f"Time {time}: The fisherman catches another goldfish, and the shark eats one.")
    print(f"There are now {goldfish_in_pond} goldfish in the pond. The basket is now full with {basket_count} fish.")
    print("The fisherman must leave for 20 minutes to get a new basket. He will return at Time 40.")
    print("---")

    # Time 30
    time = 30
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 30
    print(f"Time {time}: While the fisherman is away, the shark eats a goldfish. {goldfish_in_pond} goldfish remain.")
    print("---")
    
    # Time 40
    time = 40
    # Shark's meal
    goldfish_in_pond -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 40
    print(f"Time {time}: The fisherman returns. At the same time, the shark eats another goldfish.")
    print(f"{goldfish_in_pond} goldfish are left in the pond. The fisherman resumes fishing with Rod B.")
    print("---")
    
    # Time 45
    time = 45
    goldfish_in_pond -= 1
    # new basket has 1 fish
    print(f"Time {time}: The fisherman catches a goldfish. Only {goldfish_in_pond} goldfish is left in the pond.")
    print("---")
    
    # Time 50
    time = 50
    goldfish_in_pond -= 1
    print(f"Time {time}: The fisherman catches the last goldfish. There are now {goldfish_in_pond} free-swimming goldfish in the pond.")
    print("The shark tries to eat at this time but finds no fish.")
    print("---")
    
    # Final Calculation
    print("Now we determine when the shark dies.")
    print(f"The shark's last meal was at Time {shark_last_meal_time}.")
    print(f"The total number of fish eaten by the shark is {shark_eaten_count}.")
    print("The rule for survival extension is: 'if the shark has eaten MORE THAN 4 fish within the last 60 minutes'.")
    print("Since the shark has eaten exactly 4 fish, the extension does not apply.")
    starvation_time = 11
    print(f"The standard starvation time is {starvation_time} minutes.")
    
    final_time = shark_last_meal_time + starvation_time
    print("\nThe shark will die after 11 minutes from its last meal.")
    print(f"The final time is calculated as: {shark_last_meal_time} + {starvation_time} = {final_time} minutes.")

solve_fish_pond_problem()
<<<51>>>