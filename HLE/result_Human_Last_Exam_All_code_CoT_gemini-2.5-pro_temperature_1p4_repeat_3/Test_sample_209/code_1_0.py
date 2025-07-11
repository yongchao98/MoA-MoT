import math

def solve_fish_pond_problem():
    """
    Simulates the fish pond problem to find the earliest time it's empty.
    """
    # Initial State
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eaten_count = 0
    shark_last_meal_time = -1
    fisherman_catches_made = 0
    
    print("Finding the earliest time for the pond to be empty.")
    print("The optimal strategy involves using Rod B and letting the fisherman catch the last fish.")
    print("--- Simulation Start ---")
    print(f"Time: {time} min. Pond: {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}. Shark has eaten: {shark_eaten_count}.")

    # --- Phase 1: Fishing until the basket is full ---
    # Fisherman uses Rod B (5 min/fish). Basket capacity is 4.
    # Time to fill basket = 4 fish * 5 min/fish = 20 minutes.
    time_to_fill_basket = 20
    
    # In these 20 minutes, the shark (eats every 10 mins) eats at t=10 and t=20.
    shark_eats_in_phase1 = math.floor(time_to_fill_basket / 10)
    fisherman_catches_in_phase1 = 4
    
    time = time_to_fill_basket
    goldfish_in_pond -= (shark_eats_in_phase1 + fisherman_catches_in_phase1)
    fisherman_catches_made += fisherman_catches_in_phase1
    shark_eaten_count += shark_eats_in_phase1
    shark_last_meal_time = time # Shark's last meal in this phase is at t=20
    goldfish_in_basket = 4 

    print(f"\nTime: {time} min. Fisherman uses Rod B for 20 minutes.")
    print(f"Fisherman catches {fisherman_catches_in_phase1} fish. Shark eats {shark_eats_in_phase1} fish.")
    print(f"Pond: {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket} (Full!). Shark has eaten: {shark_eaten_count}.")
    print("Fisherman must now make a 20-minute trip for a new basket.")

    # --- Phase 2: Fisherman's basket trip ---
    trip_duration = 20
    time_of_return = time + trip_duration
    
    # Shark continues to eat every 10 mins (at t=30, t=40)
    shark_eats_during_trip = math.floor(time_of_return / 10) - math.floor(time / 10)
    
    goldfish_in_pond -= shark_eats_during_trip
    shark_eaten_count += shark_eats_during_trip
    shark_last_meal_time = time_of_return # Shark's last meal is at t=40

    time = time_of_return
    goldfish_in_basket = 0

    print(f"\nTime: {time} min. Fisherman returns from the 20-minute basket trip.")
    print(f"While he was away, the shark ate {shark_eats_during_trip} more fish.")
    print(f"Pond: {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}. Shark has eaten: {shark_eaten_count}.")
    print(f"The shark's last meal was at t={shark_last_meal_time} min.")

    # --- Phase 3: Fishing after trip until pond is empty of goldfish ---
    goldfish_remaining = goldfish_in_pond
    time_to_catch_remaining = goldfish_remaining * 5

    # Fisherman starts fishing at t=40. Catches are at t=45, t=50.
    # Shark's next meal is scheduled for t=50.
    # To ensure the earliest finish, we let the fisherman catch the last fish at t=50.
    # This prevents the shark from eating again.
    time_pond_empty = time + time_to_catch_remaining
    
    fisherman_catches_made += goldfish_remaining
    goldfish_in_pond = 0
    
    print(f"\nTime: {time_pond_empty} min. Fisherman fishes the remaining {goldfish_remaining} goldfish.")
    print("At t=50 min, the fisherman catches the last fish, preventing the shark from eating.")
    print(f"Pond: {goldfish_in_pond} goldfish. All goldfish are now gone.")
    print(f"The shark's last meal remains at t={shark_last_meal_time} min.")

    # --- Phase 4: Shark starvation ---
    # Shark survival rule: 11 minutes, unless it has eaten MORE THAN 4 fish in the last 60 mins.
    # Total eaten by shark is 4. This is not > 4, so the standard 11-minute timer applies.
    shark_survival_window = 11
    time_of_shark_death = shark_last_meal_time + shark_survival_window

    print("\n--- Final Calculation ---")
    print(f"The pond is empty of goldfish at t={time_pond_empty} minutes.")
    print(f"The shark last ate at t={shark_last_meal_time} minutes.")
    print(f"The shark has eaten a total of {shark_eaten_count} fish. This is not 'more than 4', so its survival window is 11 minutes.")
    print("The shark will starve 11 minutes after its last meal.")
    print(f"Final Time = (Time of Shark's Last Meal) + (Shark Survival Window)")
    print(f"Final Time = {shark_last_meal_time} + {shark_survival_window} = {time_of_shark_death}")
    print("\nTherefore, the earliest possible time when there are no more fish left in the pond is 51 minutes.")

solve_fish_pond_problem()