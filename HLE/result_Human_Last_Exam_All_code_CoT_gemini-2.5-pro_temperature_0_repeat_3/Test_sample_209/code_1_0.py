def solve_fish_pond_problem():
    """
    Calculates the earliest time the pond is empty by simulating the optimal strategy.
    """
    # --- Initial State ---
    pond_goldfish = 10
    basket_goldfish = 0
    shark_eaten_count = 0
    shark_last_meal_time = 0
    shark_is_powered_up = False

    print("Timeline of Events to Empty the Pond:")
    print("------------------------------------")
    print("t=0: The fisherman starts with Rod B (catches a fish every 5 minutes).")
    print(f"       State: Pond={pond_goldfish} goldfish, Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # --- Fisherman uses Rod B ---
    # t=5
    pond_goldfish -= 1
    basket_goldfish += 1
    print("t=5: Fisherman catches a goldfish with Rod B.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}\n")

    # t=10
    pond_goldfish -= 1
    basket_goldfish += 1
    # Shark eats from the pond
    pond_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 10
    print("t=10: Fisherman catches a goldfish with Rod B. The shark eats a free-swimming goldfish.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=15
    pond_goldfish -= 1
    basket_goldfish += 1
    print("t=15: Fisherman catches a goldfish with Rod B. The basket now has 3 fish.")
    print("        To avoid the basket overflowing on the next catch, he switches to Rod A and starts feeding the shark.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}\n")

    # --- Fisherman uses Rod A and feeds the shark ---
    # t=17 (Feeding takes 2 minutes)
    basket_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 17
    print("t=17: Fisherman finishes feeding the shark. This resets the shark's 10-minute eating timer.")
    print(f"       State: Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=22 (Rod A catches a fish every 7 minutes, 15+7=22)
    pond_goldfish -= 1
    basket_goldfish += 1
    print("t=22: Fisherman catches a goldfish with Rod A and immediately starts feeding the shark again.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}\n")

    # t=24
    basket_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 24
    print("t=24: Fisherman finishes feeding the shark.")
    print(f"       State: Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=29 (22+7=29)
    pond_goldfish -= 1
    basket_goldfish += 1
    print("t=29: Fisherman catches a goldfish with Rod A and starts feeding the shark.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}\n")

    # t=31
    basket_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 31
    print("t=31: Fisherman finishes feeding the shark.")
    print(f"       State: Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=34 (Shark's natural meal, 10 mins after its last meal at t=24, before being fed at t=31)
    pond_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 34
    print("t=34: The shark eats a free-swimming goldfish on its own.")
    print(f"       State: Pond={pond_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=36 (29+7=36)
    pond_goldfish -= 1
    basket_goldfish += 1
    print("t=36: Fisherman catches a goldfish with Rod A and starts the final feeding.")
    print(f"       State: Pond={pond_goldfish}, Basket={basket_goldfish}\n")

    # t=38
    basket_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 38
    shark_is_powered_up = True
    print("t=38: Fisherman finishes feeding the shark. This is the shark's 6th meal!")
    print("        The shark is now powered up and will eat every 2 minutes.")
    print(f"       State: Basket={basket_goldfish}, Shark Meals={shark_eaten_count}\n")

    # --- Powered-up shark clears the pond ---
    # t=40 (38+2=40)
    pond_goldfish -= 1
    shark_eaten_count += 1
    shark_last_meal_time = 40
    print("t=40: The powered-up shark eats a goldfish.")
    print(f"       State: Pond={pond_goldfish}, Shark Meals={shark_eaten_count}\n")

    # t=42 (40+2=42)
    pond_goldfish -= 1
    shark_eaten_count += 1
    time_pond_empty = 42
    shark_last_meal_time = time_pond_empty
    print(f"t={time_pond_empty}: The shark eats the last goldfish from the pond.")
    print(f"       State: Pond={pond_goldfish}, Shark Meals={shark_eaten_count}. All goldfish are gone!\n")

    # --- Final Calculation ---
    print("------------------------------------")
    print("Final Calculation:")
    print(f"The pond is empty of goldfish at t={time_pond_empty} minutes.")
    
    shark_starvation_base = 11
    shark_starvation_bonus = 4
    # The shark has eaten > 4 fish in the last 60 minutes, so its survival time is extended.
    total_starvation_time = shark_starvation_base + shark_starvation_bonus
    
    print(f"The shark's last meal was at t={shark_last_meal_time}. Because it has eaten more than 4 fish recently,")
    print(f"its survival time is extended by {shark_starvation_bonus} minutes, from {shark_starvation_base} to {total_starvation_time} minutes.")
    
    final_time = time_pond_empty + total_starvation_time
    
    print("\nThe shark starves and disappears. The pond is now completely empty of fish.")
    print("\nFinal Equation:")
    print(f"{time_pond_empty} + ({shark_starvation_base} + {shark_starvation_bonus}) = {final_time}")


solve_fish_pond_problem()
<<<57>>>