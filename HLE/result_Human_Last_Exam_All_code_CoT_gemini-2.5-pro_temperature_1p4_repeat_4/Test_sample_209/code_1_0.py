import math

def solve():
    """
    Simulates the fisherman and shark problem to find the earliest time all fish are gone.
    """
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    basket_capacity = 4

    shark_meals_total = 0
    shark_meal_timestamps = []
    shark_eat_rate = 10  # minutes

    # Fisherman state
    # Rod A: 7 min/catch, 1-handed
    # Rod B: 5 min/catch, 2-handed
    fisherman_rod = 'B'
    time_since_rod_A_started = -1
    rod_b_continuous_use_time = 0

    # Event timers
    next_fisher_catch = 5
    next_shark_eat = 10
    
    # Track events for logging
    event_log = []

    # Main simulation loop
    while goldfish_in_pond > 0:
        time_to_next_event = min(next_fisher_catch, next_shark_eat)
        time += time_to_next_event

        # Update timers
        next_fisher_catch -= time_to_next_event
        next_shark_eat -= time_to_next_event
        if fisherman_rod == 'B':
            rod_b_continuous_use_time += time_to_next_event
        elif fisherman_rod == 'A' and time_since_rod_A_started != -1:
             time_since_rod_A_started += time_to_next_event


        # --- Event Handling ---

        # Shark eats event
        if next_shark_eat == 0:
            if goldfish_in_pond > 0:
                goldfish_in_pond -= 1
                shark_meals_total += 1
                shark_meal_timestamps.append(time)
                event_log.append(f"Time {time}: Shark eats a free-swimming goldfish. {goldfish_in_pond} left in pond. Total shark meals: {shark_meals_total}.")

                if shark_meals_total >= 6 and shark_eat_rate == 10:
                    shark_eat_rate = 2
                    event_log.append(f"Time {time}: Shark has eaten 6 fish and becomes stronger. Now eats every 2 minutes.")
            
            next_shark_eat = shark_eat_rate

        # Fisherman event (catch or feed)
        if next_fisher_catch == 0:
            # Fisherman catches a fish
            if goldfish_in_pond > 0 and goldfish_in_basket < basket_capacity:
                goldfish_in_pond -= 1
                goldfish_in_basket += 1
                event_log.append(f"Time {time}: Fisherman catches a goldfish with Rod {fisherman_rod}. {goldfish_in_pond} left in pond, {goldfish_in_basket} in basket.")

            # Logic after a catch, or if basket was already full
            
            # --- Fisherman's Decision Logic ---
            
            # 1. Basket is full, must feed
            if goldfish_in_basket == basket_capacity:
                # Switch to Rod A if not already using it
                if fisherman_rod == 'B':
                    event_log.append(f"Time {time}: Basket is full. Fisherman must switch to Rod A to feed the shark.")
                    fisherman_rod = 'A'
                    time_since_rod_A_started = 0
                
                # Feed the shark
                feed_time = 2
                time += feed_time
                if fisherman_rod == 'A' and time_since_rod_A_started != -1:
                    time_since_rod_A_started += feed_time
                next_shark_eat -= feed_time # Shark timer continues

                goldfish_in_basket -= 1
                shark_meals_total += 1
                shark_meal_timestamps.append(time)
                next_shark_eat = shark_eat_rate # Feeding resets shark's hunger
                event_log.append(f"Time {time}: Fisherman feeds a goldfish to the shark. Basket has {goldfish_in_basket} fish. Total shark meals: {shark_meals_total}.")

                if shark_meals_total >= 6 and shark_eat_rate == 10:
                    shark_eat_rate = 2
                    event_log.append(f"Time {time}: Shark has eaten 6 fish and becomes stronger. Now eats every 2 minutes.")


            # 2. Decide on next rod
            if fisherman_rod == 'A' and time_since_rod_A_started >= 7 and goldfish_in_basket < basket_capacity:
                 event_log.append(f"Time {time}: Fisherman has used Rod A long enough and switches back to Rod B.")
                 fisherman_rod = 'B'
                 rod_b_continuous_use_time = 0
                 time_since_rod_A_started = -1


            # 3. Set next catch time
            if fisherman_rod == 'A':
                 next_fisher_catch = 7
            else: # Rod B
                 next_fisher_catch = 5

    # --- Final Calculation ---
    last_meal_time = shark_meal_timestamps[-1]
    
    meals_in_last_60_mins = 0
    for t in shark_meal_timestamps:
        if last_meal_time - 60 < t <= last_meal_time:
            meals_in_last_60_mins += 1

    base_starvation_time = 11
    extension = 4 if meals_in_last_60_mins > 4 else 0
    total_starvation_time = base_starvation_time + extension
    final_time = last_meal_time + total_starvation_time
    
    # Print the step-by-step reasoning for the final answer
    print("Finding the earliest time for all fish to be gone:")
    print("-" * 50)
    print(f"1. Time when all free-swimming goldfish are gone: {last_meal_time} minutes.")
    print(f"2. The shark's last meal was at time: {last_meal_time} minutes.")
    print(f"3. Checking shark's meals in the 60 minutes prior to its last meal...")
    print(f"   - Number of fish eaten by shark in the last 60 minutes: {meals_in_last_60_mins}.")
    if extension > 0:
        print(f"   - Since this is more than 4, the shark's survival time is extended by {extension} minutes.")
    else:
        print(f"   - Since this is not more than 4, the shark's survival time is not extended.")
    print(f"4. Calculating total starvation time:")
    print(f"   - Base starvation time: {base_starvation_time} minutes")
    print(f"   - Time extension: {extension} minutes")
    print(f"   - Final equation for shark survival: {base_starvation_time} + {extension} = {total_starvation_time} minutes.")
    print(f"5. Calculating the final time:")
    print(f"   - Final equation: {last_meal_time} (time of last meal) + {total_starvation_time} (starvation time) = {final_time} minutes.")
    print("-" * 50)
    print(f"The earliest possible time when there are no more fish left in the pond is {final_time} minutes.")

solve()
<<<60>>>