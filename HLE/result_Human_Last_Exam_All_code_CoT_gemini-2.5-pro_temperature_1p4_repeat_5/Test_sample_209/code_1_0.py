def solve_fish_pond_problem():
    """
    Simulates the optimal strategy to empty the fish pond at the earliest possible time.
    """
    # Initial state of the simulation
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eats_count = 0
    # The shark's eating schedule is every 10 minutes.
    # We can model its last meal as being at T=-10 to have the first meal happen at T=10.
    shark_last_eaten_time = -10

    # Print initial state
    print(f"Time: {time} minutes")
    print(f"The pond has {goldfish_in_pond} goldfish and 1 shark.")
    print("The fisherman's basket is empty (0/4).")
    print("Optimal strategy: Use Rod B as much as possible.")
    print("---")

    # --- Event 1: T=5 min ---
    time = 5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time: {time} minutes")
    print("Fisherman catches a goldfish with Rod B (5-minute action).")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print("---")

    # --- Event 2: T=10 min ---
    time = 10
    # Fisherman's action
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time: {time} minutes")
    print("Fisherman catches another goldfish with Rod B.")
    # Shark's action
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_eaten_time = time
    print("Simultaneously, the shark eats a free-swimming goldfish.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print(f"Shark has eaten {shark_eats_count} goldfish so far.")
    print("---")

    # --- Event 3: T=15 min ---
    time = 15
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time: {time} minutes")
    print("Fisherman catches another goldfish with Rod B.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print("---")
    
    # --- Event 4: T=20 min (Basket becomes full) ---
    time = 20
    # Fisherman's action
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time: {time} minutes")
    print("Fisherman catches another goldfish, filling the basket.")
    # Shark's action
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_eaten_time = time
    print("Simultaneously, the shark eats another free-swimming goldfish.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print("The fisherman must leave for 20 minutes to get a new basket.")
    print("---")

    # --- Fisherman is away (T=20 to T=40) ---
    
    # --- Event 5: T=30 min ---
    time = 30
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_eaten_time = time
    print(f"Time: {time} minutes")
    print("While the fisherman is away, the shark eats another goldfish.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Shark has eaten {shark_eats_count} goldfish.")
    print("---")
    
    # --- Event 6: T=40 min (Fisherman returns) ---
    time = 40
    # Shark's action
    goldfish_in_pond -= 1
    shark_eats_count += 1
    shark_last_eaten_time = time
    print(f"Time: {time} minutes")
    print("The shark eats another goldfish.")
    # Fisherman's action
    goldfish_in_basket = 0
    print("The fisherman returns with a new empty basket and resumes fishing.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print(f"Shark has eaten {shark_eats_count} goldfish so far.")
    print("---")

    # --- Event 7: T=45 min ---
    time = 45
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time: {time} minutes")
    print("Fisherman catches another goldfish with Rod B.")
    print(f"Pond now has {goldfish_in_pond} goldfish. Basket: {goldfish_in_basket}/4.")
    print("---")

    # --- Event 8: T=50 min (Final goldfish) ---
    time = 50
    print(f"Time: {time} minutes")
    print("One goldfish remains. To get the earliest time, the fisherman must catch it before the shark does.")
    goldfish_in_pond -= 1
    print("Fisherman catches the last goldfish.")
    print(f"All {goldfish_in_pond} goldfish are now gone from the pond.")
    print("---")
    
    # --- Final Calculation: Shark's fate ---
    print("Calculating the final time when the pond is completely empty.")
    # The rule is: if shark ate > 4 fish in the last 60 minutes, survival is 11+4=15 mins.
    # At T=50, the shark's meals were at T=10, 20, 30, 40. That's exactly 4 meals.
    # The condition "> 4" is false.
    shark_survival_time = 11
    print(f"The shark has eaten {shark_eats_count} fish, which is not more than 4 in the last 60 minutes.")
    print(f"Therefore, its survival time without food is the standard {shark_survival_time} minutes.")
    print(f"The shark's last meal was at T={shark_last_eaten_time} minutes.")
    
    final_time = shark_last_eaten_time + shark_survival_time
    print("\nThe pond is empty of all fish once the shark starves.")
    print(f"Time of last fish disappearing = Shark's Last Meal Time + Shark's Survival Time")
    print(f"Earliest possible time = {shark_last_eaten_time} + {shark_survival_time} = {final_time}")
    
    return final_time

if __name__ == "__main__":
    solve_fish_pond_problem()
    final_answer = 51
    print(f"\n<<<{final_answer}>>>")