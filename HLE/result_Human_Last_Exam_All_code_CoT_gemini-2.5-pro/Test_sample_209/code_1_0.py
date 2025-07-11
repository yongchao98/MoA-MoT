def solve_fish_pond_problem():
    """
    This function calculates the earliest time the pond is empty by simulating the optimal strategy.
    The optimal strategy is for the fisherman to use the faster Rod B, even though this involves
    a 20-minute trip to get a new basket. This clears the goldfish faster than any other method.
    """

    # Initial state
    initial_goldfish = 10
    fisherman_catches = 0
    shark_eats = 0

    print("--- Finding the Optimal Timeline ---")
    print(f"Initial state: {initial_goldfish} goldfish, 1 shark.\n")

    # Phase 1: Fisherman uses Rod B until the basket is full.
    print("Phase 1 (Time 0-20 mins): Fisherman uses Rod B (1 fish/5 mins), Shark eats (1 fish/10 mins)")
    # Fish removed by t=20: 4 by fisherman (at t=5, 10, 15, 20) and 2 by shark (at t=10, 20).
    fisherman_catches_p1 = 4
    shark_eats_p1 = 2
    goldfish_removed_p1 = fisherman_catches_p1 + shark_eats_p1
    goldfish_remaining_p1 = initial_goldfish - goldfish_removed_p1
    
    print(f"At t=5, fisherman catches 1 fish. (9 left)")
    print(f"At t=10, fisherman catches 1 fish and shark eats 1. (7 left)")
    print(f"At t=15, fisherman catches 1 fish. (6 left)")
    print(f"At t=20, fisherman catches his 4th fish and shark eats its 2nd. The basket is full.")
    print(f"Total fish removed by t=20: {fisherman_catches_p1} (caught) + {shark_eats_p1} (eaten) = {goldfish_removed_p1}")
    print(f"Goldfish remaining at t=20: {goldfish_remaining_p1}\n")
    
    # Phase 2: Fisherman is on a 20-minute basket trip.
    print("Phase 2 (Time 20-40 mins): Fisherman on 20-min basket trip")
    # Fish removed during trip: 2 by shark (at t=30, 40).
    shark_eats_p2 = 2
    goldfish_removed_p2 = shark_eats_p2
    goldfish_remaining_p2 = goldfish_remaining_p1 - goldfish_removed_p2
    shark_last_meal_time = 40
    
    print(f"At t=30, the shark eats a goldfish. (3 left)")
    print(f"At t=40, the shark eats another goldfish. (2 left)")
    print(f"The fisherman returns at t=40. The shark's last meal was at t={shark_last_meal_time}.\n")

    # Phase 3: Fisherman clears the remaining fish.
    print("Phase 3 (Time 40-50 mins): Fisherman returns and catches the remaining fish.")
    # Fish removed: 2 by fisherman (at t=45, 50).
    fisherman_catches_p3 = 2
    goldfish_remaining_p3 = goldfish_remaining_p2 - fisherman_catches_p3
    
    print(f"At t=45, the fisherman catches a fish. (1 left)")
    print(f"At t=50, the fisherman catches the last goldfish. ({goldfish_remaining_p3} left)")
    print("All goldfish are gone from the pond at t=50.\n")

    # Phase 4: Calculate when the shark dies.
    print("Phase 4: Calculate the shark's death.")
    # Total shark meals: 2 in Phase 1 + 2 in Phase 2 = 4 meals.
    shark_total_meals = shark_eats_p1 + shark_eats_p2
    # The condition for extended survival is "> 4 meals in the last 60 mins".
    # At t=50, the shark has eaten 4 meals in the last 60 mins. 4 is not > 4.
    survival_time_extension_threshold = 4
    shark_survival_time = 11
    
    print(f"The shark's last meal was at time t = {shark_last_meal_time}.")
    print(f"The shark ate a total of {shark_total_meals} fish.")
    print(f"The survival time extension applies if the shark eats more than {survival_time_extension_threshold} fish in 60 minutes.")
    print(f"Since {shark_total_meals} is not greater than {survival_time_extension_threshold}, the standard survival time of {shark_survival_time} minutes applies.")
    
    # Final Equation
    shark_death_time = shark_last_meal_time + shark_survival_time
    print("\n--- Final Calculation ---")
    print("The final event is the death of the shark, which happens after all goldfish are gone.")
    print("Shark Death Time = Time of Last Meal + Survival Duration")
    print(f"Final equation: {shark_last_meal_time} + {shark_survival_time} = {shark_death_time}")

    final_answer = shark_death_time
    print(f"\nThe earliest possible time when there are no more fish left in the pond is {final_answer} minutes.")

solve_fish_pond_problem()