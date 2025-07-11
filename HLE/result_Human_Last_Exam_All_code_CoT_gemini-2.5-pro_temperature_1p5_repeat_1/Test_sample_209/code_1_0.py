import math

def solve_fish_problem():
    """
    Simulates the fish pond problem to find the earliest time it becomes empty.

    The optimal strategy simulated here is for the fisherman to use his fastest
    rod (Rod B) continuously. This fills the basket quickly, leading to a
    20-minute trip for a new basket, but this proves to be one of the
    fastest ways to clear the pond.
    """
    # Initial state
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    shark_eats_counter = 0
    fisherman_catches_counter = 0

    # Timers for next events
    fisherman_next_catch_time = 5  # Using Rod B (5 mins)
    shark_next_eat_time = 10       # Shark's initial interval

    # Fisherman state
    rod_b_time_used = 0
    on_basket_trip_until = 0

    # Log of events to be printed
    events = []

    def log_event(t, description):
        events.append(f"Time {t:2}: {description}")

    log_event(time, f"Start. Pond: {goldfish_in_pond} goldfish. Fisherman starts using Rod B.")

    while goldfish_in_pond > 0:
        # Check if the fisherman's basket trip is over
        if on_basket_trip_until > 0 and time >= on_basket_trip_until:
            log_event(time, "Fisherman returns with a new basket. Fishing can resume.")
            goldfish_in_basket = 0
            on_basket_trip_until = 0
            # Schedule next catch now that he is back
            fisherman_next_catch_time = time + 5

        # Determine the time of the very next event
        if on_basket_trip_until > 0: # Fisherman is away, only shark can act
             next_event_time = shark_next_eat_time
        else:
            next_event_time = min(fisherman_next_catch_time, shark_next_eat_time)
            if next_event_time == float('inf'):
                break

        # Advance time to the next event
        time = next_event_time
        
        if time > 200: # Safety break to prevent infinite loops
            log_event(time, "Simulation timeout.")
            break

        # Process events happening at the current time
        # Note: Both events can happen at the same time (e.g., at t=10)
        
        # Check for Fisherman's Catch Event
        if time == fisherman_next_catch_time and goldfish_in_pond > 0 and on_basket_trip_until == 0:
            goldfish_in_pond -= 1
            goldfish_in_basket += 1
            fisherman_catches_counter += 1
            rod_b_time_used += 5
            log_event(time, f"Fisherman catches a goldfish. Pond: {goldfish_in_pond}, Basket: {goldfish_in_basket}.")
            
            # Schedule the next catch if the pond is not empty
            if goldfish_in_pond > 0:
                fisherman_next_catch_time = time + 5

            # Check if the basket is full
            if goldfish_in_basket == 4:
                log_event(time, "Basket is full! Fisherman leaves for 20 minutes to get a new one.")
                on_basket_trip_until = time + 20
                fisherman_next_catch_time = float('inf') # Cannot fish while away

        # Check for Shark's Eating Event
        if time == shark_next_eat_time and goldfish_in_pond > 0:
            goldfish_in_pond -= 1
            shark_eats_counter += 1
            log_event(time, f"Shark eats a goldfish. Pond: {goldfish_in_pond}, Shark has eaten {shark_eats_counter}.")
            
            # Schedule the shark's next meal
            shark_next_eat_time = time + 10 # Shark is not powered-up in this scenario

    # Final reporting
    final_time = time
    
    print("--- Simulation Timeline ---")
    for event in events:
        print(event)
    print("--- End of Simulation ---\n")
    
    print("--- Final Calculation ---")
    print("The final time is determined by the moment the last fish is removed from the pond.")
    print("Below is the step-by-step removal of all 10 goldfish:")
    print("Initial state at t=0: 10 goldfish")
    print("1.  t=5:  Fisherman catches 1. Pond: 10 - 1 = 9")
    print("2.  t=10: Fisherman catches 1. Pond: 9 - 1 = 8")
    print("3.  t=10: Shark eats 1.        Pond: 8 - 1 = 7")
    print("4.  t=15: Fisherman catches 1. Pond: 7 - 1 = 6")
    print("5.  t=20: Fisherman catches 1. Pond: 6 - 1 = 5 (Basket is now full)")
    print("6.  t=20: Shark eats 1.        Pond: 5 - 1 = 4")
    print("7.  t=30: Shark eats 1.        Pond: 4 - 1 = 3 (Fisherman is away)")
    print("8.  t=40: Shark eats 1.        Pond: 3 - 1 = 2 (Fisherman returns)")
    print("9.  t=45: Fisherman catches 1. Pond: 2 - 1 = 1")
    print("10. t=50: Shark eats 1.        Pond: 1 - 1 = 0")
    print(f"\nThe pond is empty at the time of the last event. The earliest possible time is {final_time} minutes.")

solve_fish_problem()