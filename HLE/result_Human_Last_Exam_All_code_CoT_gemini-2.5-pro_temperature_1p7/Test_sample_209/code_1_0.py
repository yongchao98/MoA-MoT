import sys

def solve():
    """
    This function calculates the earliest time the pond is empty of all fish
    by simulating the optimal strategy for the fisherman.
    """

    # Initial state of the pond and fisherman
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eaten_total = 0
    time = 0

    print("--- Solving the Fish Pond Puzzle ---")
    print("The goal is to find the earliest time when no fish (goldfish or shark) are left in the pond.")
    print("The fisherman's best strategy is to use the fastest rod (Rod B) to clear the goldfish as quickly as possible.")
    print("\n--- Simulation Start ---")
    print(f"Initial State: Time = {time} mins, Goldfish in Pond = {goldfish_in_pond}\n")

    # Time: 5 minutes
    # The fisherman uses Rod B, which catches a fish every 5 minutes.
    time = 5
    fisherman_catch_1 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= fisherman_catch_1
    goldfish_in_basket += fisherman_catch_1
    print(f"Time = {time} mins: The fisherman catches his first goldfish using Rod B.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {fisherman_catch_1} (caught) = {goldfish_in_pond} goldfish remaining.")
    print(f"Fisherman's basket: {goldfish_in_basket} goldfish.\n")

    # Time: 10 minutes
    # The shark eats a fish (every 10 mins), and the fisherman catches another.
    time = 10
    shark_eat_1 = 1
    fisherman_catch_2 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= (shark_eat_1 + fisherman_catch_2)
    goldfish_in_basket += fisherman_catch_2
    shark_eaten_total += shark_eat_1
    print(f"Time = {time} mins: The shark eats a goldfish and the fisherman simultaneously catches another.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {shark_eat_1} (shark) - {fisherman_catch_2} (caught) = {goldfish_in_pond} goldfish remaining.")
    print(f"Fisherman's basket: {goldfish_in_basket} goldfish.\n")

    # Time: 15 minutes
    time = 15
    fisherman_catch_3 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= fisherman_catch_3
    goldfish_in_basket += fisherman_catch_3
    print(f"Time = {time} mins: The fisherman catches his third goldfish.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {fisherman_catch_3} (caught) = {goldfish_in_pond} goldfish remaining.")
    print(f"Fisherman's basket: {goldfish_in_basket} goldfish.\n")

    # Time: 20 minutes
    # The basket becomes full, forcing the fisherman to take action.
    time = 20
    shark_eat_2 = 1
    fisherman_catch_4 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= (shark_eat_2 + fisherman_catch_4)
    goldfish_in_basket += fisherman_catch_4
    shark_eaten_total += shark_eat_2
    print(f"Time = {time} mins: The shark eats another, and the fisherman catches his fourth.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {shark_eat_2} (shark) - {fisherman_catch_4} (caught) = {goldfish_in_pond} goldfish remaining.")
    print(f"Fisherman's basket is now full with {goldfish_in_basket} goldfish.\n")

    # The 20-minute trip for a new basket
    basket_trip_time = 20
    time_fisherman_leaves = time
    time_fisherman_returns = time_fisherman_leaves + basket_trip_time
    print(f"--- Fisherman's Trip for New Basket ---")
    print(f"At T={time_fisherman_leaves} mins, the fisherman must get a new basket, which takes {basket_trip_time} minutes.")
    print(f"Equation: {time_fisherman_leaves} (leaves) + {basket_trip_time} (trip) = {time_fisherman_returns} (returns)")
    
    # Events during the fisherman's absence
    time = 30
    shark_eat_3 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= shark_eat_3
    shark_eaten_total += shark_eat_3
    print(f"Time = {time} mins: While the fisherman is away, the shark eats.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {shark_eat_3} (shark) = {goldfish_in_pond} goldfish remaining.")

    time = 40
    shark_eat_4 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= shark_eat_4
    shark_eaten_total += shark_eat_4
    time_of_sharks_last_meal = time
    print(f"Time = {time} mins: The shark eats again just as the fisherman returns with a new, empty basket.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {shark_eat_4} (shark) = {goldfish_in_pond} goldfish remaining.")
    print(f"The shark has now eaten {shark_eaten_total} goldfish. Its last meal was at T={time_of_sharks_last_meal}.\n")
    goldfish_in_basket = 0 # Basket is reset

    # Fisherman clears the remaining goldfish
    print(f"--- Clearing the Last Goldfish ---")
    time = 45
    fisherman_catch_5 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= fisherman_catch_5
    goldfish_in_basket += fisherman_catch_5
    print(f"Time = {time} mins: The fisherman resumes fishing and catches another goldfish.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {fisherman_catch_5} (caught) = {goldfish_in_pond} goldfish remaining.\n")

    time = 50
    fisherman_catch_6 = 1
    goldfish_in_pond_before = goldfish_in_pond
    goldfish_in_pond -= fisherman_catch_6
    print(f"Time = {time} mins: The fisherman catches the final goldfish.")
    print(f"Equation: {goldfish_in_pond_before} (pond) - {fisherman_catch_6} (caught) = {goldfish_in_pond} goldfish remaining.")
    print("All goldfish have now been removed from the pond.\n")

    # Final calculation for the shark's demise
    shark_starvation_time = 11
    print(f"--- Final Calculation ---")
    print("The final step is for the shark to perish. This happens when it doesn't eat for a specific duration.")
    print(f"The shark's last meal was at T={time_of_sharks_last_meal} mins. The base starvation time is {shark_starvation_time} minutes.")
    print(f"The survival extension (for eating > 4 fish in 60 mins) is not triggered because the shark ate exactly {shark_eaten_total} fish, not more than 4.")
    
    final_time = time_of_sharks_last_meal + shark_starvation_time
    print("\nFinal Time = Time of Shark's Last Meal + Shark Starvation Time")
    print(f"Final Time = {time_of_sharks_last_meal} + {shark_starvation_time} = {final_time}")
    
    return final_time

if __name__ == '__main__':
    # The result is printed by the function itself.
    final_answer = solve()
    # To conform to the output format, we also print the final answer bracketed at the end.
    sys.stdout.flush() # ensure all prints are out
    print(f"\n<<<{final_answer}>>>")