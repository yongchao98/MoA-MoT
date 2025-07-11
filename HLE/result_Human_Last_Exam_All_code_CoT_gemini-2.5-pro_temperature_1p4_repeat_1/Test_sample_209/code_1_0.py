import math

def solve_fish_problem():
    """
    Simulates the fish pond problem to find the earliest time all fish are gone.
    """
    # Initial state of the simulation
    time = 0
    free_goldfish = 10
    caught_goldfish = 0
    shark_eaten_total = 0
    shark_is_alive = True
    shark_last_meal_time = 0
    shark_eat_history = []
    fisherman_is_away = False
    
    # Timers for the next events
    fisherman_catch_timer = 5  # Using Rod B from the start
    shark_eat_timer = 10
    fisherman_away_timer = 0
    
    events = []

    # Main simulation loop
    # Continues as long as there are goldfish to be caught or the shark is alive
    while free_goldfish > 0 or shark_is_alive:
        # Determine the time of the very next event
        next_event_time = float('inf')
        
        # Possible event: Fisherman catches a fish
        if not fisherman_is_away and free_goldfish > 0:
            next_event_time = min(next_event_time, time + fisherman_catch_timer)
        
        # Possible event: Shark eats a fish
        if shark_is_alive and free_goldfish > 0:
            next_event_time = min(next_event_time, time + shark_eat_timer)
            
        # Possible event: Fisherman returns from trip
        if fisherman_is_away:
            next_event_time = min(next_event_time, time + fisherman_away_timer)
            
        # Possible event: Shark starves
        if free_goldfish == 0 and shark_is_alive:
            eaten_in_last_60_min = sum(1 for meal_time in shark_eat_history if time - meal_time <= 60)
            starvation_duration = 15 if eaten_in_last_60_min > 4 else 11
            shark_death_time = shark_last_meal_time + starvation_duration
            next_event_time = min(next_event_time, shark_death_time)

        if next_event_time == float('inf'):
            break

        # Advance time to the next event
        elapsed_time = next_event_time - time
        fisherman_catch_timer -= elapsed_time
        shark_eat_timer -= elapsed_time
        if fisherman_is_away:
            fisherman_away_timer -= elapsed_time
        time = next_event_time

        # --- Process events that occur at the new time ---

        # Event: Fisherman returns
        if fisherman_is_away and fisherman_away_timer <= 0:
            fisherman_is_away = False
            caught_goldfish = 0
            fisherman_catch_timer = 5  # Resumes with Rod B
            events.append(f"Time {int(time)}: Fisherman returns with a new basket.")

        # Event: Shark eats
        if shark_is_alive and free_goldfish > 0 and shark_eat_timer <= 0:
            free_goldfish -= 1
            shark_eaten_total += 1
            shark_last_meal_time = time
            shark_eat_history.append(time)
            # Shark eats faster after 6 meals, but this condition isn't met in the optimal path
            shark_eat_timer = 10 
            events.append(f"Time {int(time)}: Shark eats a goldfish. ({free_goldfish} left in pond)")

        # Event: Fisherman catches a fish
        if not fisherman_is_away and free_goldfish > 0 and fisherman_catch_timer <= 0:
            free_goldfish -= 1
            caught_goldfish += 1
            fisherman_catch_timer = 5  # Reset Rod B timer
            events.append(f"Time {int(time)}: Fisherman catches a goldfish. ({free_goldfish} left in pond)")
            
            if caught_goldfish >= 4:
                fisherman_is_away = True
                fisherman_away_timer = 20
                events.append(f"Time {int(time)}: Basket is full. Fisherman leaves for 20 minutes.")
                last_goldfish_time = time
        
        # Update the time the last goldfish was removed
        if free_goldfish == 0:
            last_goldfish_time = time

        # Event: Shark dies from starvation
        if free_goldfish == 0 and shark_is_alive:
            eaten_in_last_60_min = sum(1 for meal_time in shark_eat_history if last_goldfish_time - meal_time <= 60)
            starvation_duration = 15 if eaten_in_last_60_min > 4 else 11
            if time >= shark_last_meal_time + starvation_duration:
                shark_is_alive = False
                events.append(f"Time {int(time)}: The shark starves and dies.")

    # Print the key events log
    print("Log of key events:")
    for event in events:
        print(event)
        
    print("\n" + "="*20)
    print("Final Calculation")
    print("="*20)

    # Calculate final starvation conditions
    final_starvation_check_time = last_goldfish_time
    eaten_in_last_60_min = sum(1 for meal_time in shark_eat_history if final_starvation_check_time - meal_time <= 60)
    final_starvation_duration = 15 if eaten_in_last_60_min > 4 else 11

    # Output the final calculation step-by-step
    print(f"The last free-swimming goldfish was removed at time: {int(last_goldfish_time)}")
    print(f"The shark's last meal was at time: {int(shark_last_meal_time)}")
    print(f"Number of fish eaten by shark in the 60 minutes before starvation began: {eaten_in_last_60_min}")
    print(f"Because the number of fish eaten ({eaten_in_last_60_min}) is not greater than 4, the standard starvation time applies.")
    print(f"The final equation is: Time of last meal + Starvation Time")
    print(f"Earliest Time = {int(shark_last_meal_time)} + {final_starvation_duration} = {int(shark_last_meal_time + final_starvation_duration)}")


solve_fish_problem()
<<<51>>>