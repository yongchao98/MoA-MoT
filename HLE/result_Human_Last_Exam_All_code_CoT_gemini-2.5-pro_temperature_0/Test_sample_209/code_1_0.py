def solve_fish_pond_problem():
    """
    Simulates the optimal strategy to clear the pond of all fish.
    """
    # Initial State
    goldfish_in_pond = 10
    fisherman_basket = 0
    shark_meals = 0
    shark_last_meal_time = 0
    
    print("Timeline of Events:")
    print("T=0: Fisherman starts fishing with Rod B (catches a fish every 5 minutes).")
    print(f"       Initial state: {goldfish_in_pond} goldfish in pond.\n")

    # Phase 1: Fisherman fishes until basket is full
    # T=5
    time = 5
    fisherman_basket += 1
    goldfish_in_pond -= 1
    print(f"T={time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.")

    # T=10
    time = 10
    fisherman_basket += 1
    goldfish_in_pond -= 1
    print(f"T={time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.")
    # Shark eats at the same time
    shark_meals += 1
    goldfish_in_pond -= 1
    shark_last_meal_time = 10
    print(f"T={time}: Shark eats a free-swimming goldfish. Pond: {goldfish_in_pond}, Shark Meals: {shark_meals}.\n")

    # T=15
    time = 15
    fisherman_basket += 1
    goldfish_in_pond -= 1
    print(f"T={time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.\n")

    # T=20
    time = 20
    fisherman_basket += 1
    goldfish_in_pond -= 1
    print(f"T={time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.")
    # Shark eats at the same time
    shark_meals += 1
    goldfish_in_pond -= 1
    shark_last_meal_time = 20
    print(f"T={time}: Shark eats a free-swimming goldfish. Pond: {goldfish_in_pond}, Shark Meals: {shark_meals}.")
    print(f"T={time}: Fisherman's basket is full (4 fish).\n")

    # Phase 2: Fisherman's 20-minute trip for a new basket
    print(f"T=20 to T=40: Fisherman takes a 20-minute trip to get a new basket.")
    
    # T=30
    time = 30
    shark_meals += 1
    goldfish_in_pond -= 1
    shark_last_meal_time = 30
    print(f"T={time}: While fisherman is away, shark eats. Pond: {goldfish_in_pond}, Shark Meals: {shark_meals}.")

    # T=40
    time = 40
    shark_meals += 1
    goldfish_in_pond -= 1
    shark_last_meal_time = 40
    print(f"T={time}: While fisherman is away, shark eats. Pond: {goldfish_in_pond}, Shark Meals: {shark_meals}.")
    fisherman_basket = 0
    print(f"T={time}: Fisherman returns with a new empty basket.\n")

    # Phase 3: Fisherman clears the remaining goldfish
    # T=45
    time = 45
    fisherman_basket += 1
    goldfish_in_pond -= 1
    print(f"T={time}: Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.")

    # T=50
    time = 50
    fisherman_basket += 1
    goldfish_in_pond -= 1
    t_goldfish_gone = 50
    print(f"T={time}: Fisherman catches the last goldfish. Pond: {goldfish_in_pond}, Basket: {fisherman_basket}.")
    print(f"T={t_goldfish_gone}: All 10 goldfish are gone from the pond.\n")

    # Phase 4: Calculate the final time when the shark is also gone
    print("Final Calculation:")
    print(f"The shark's last meal was at T={shark_last_meal_time}.")
    print(f"Total meals for the shark: {shark_meals}.")

    starvation_time = 11
    # Check if starvation time is extended
    if shark_meals > 4:
        starvation_time += 4
        print(f"Since the shark ate more than 4 fish, its survival time is extended by 4 minutes to {starvation_time} minutes.")
    else:
        print(f"Since the shark has not eaten more than 4 fish ({shark_meals}), its survival time is the standard {starvation_time} minutes.")

    shark_death_time = shark_last_meal_time + starvation_time
    print(f"The shark will die from starvation at T = {shark_last_meal_time} + {starvation_time} = {shark_death_time}.")

    final_time = max(t_goldfish_gone, shark_death_time)
    print(f"\nThe earliest time when no fish are left is the maximum of when the goldfish are gone (T={t_goldfish_gone}) and when the shark is gone (T={shark_death_time}).")
    print(f"Final Time = max({t_goldfish_gone}, {shark_death_time}) = {final_time} minutes.")

    return final_time

if __name__ == '__main__':
    final_answer = solve_fish_pond_problem()
    print(f"\n<<< {final_answer} >>>")
