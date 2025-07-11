import math

def solve_fish_puzzle():
    """
    Simulates the fish pond problem to find the earliest time all fish are gone.
    The strategy is to make the shark stronger as quickly as possible.
    """
    # --- State Variables ---
    time = 0
    goldfish_in_pond = 10
    goldfish_in_basket = 0
    
    shark_eaten_count = 0
    shark_is_strong = False
    shark_last_eaten_time = 0
    shark_eat_times = [] # a log of timestamps when the shark ate

    # --- Event Simulation ---
    print("Starting the simulation to find the earliest time the pond is empty.")
    print("-" * 60)

    # Phase 1: Fisherman catches the first fish with the faster Rod B.
    time = 5
    goldfish_in_pond -= 1
    goldfish_in_basket += 1
    print(f"Time {time:2}: Fisherman uses Rod B and catches a goldfish. {goldfish_in_pond} goldfish left in pond, {goldfish_in_basket} in basket.")
    print("          Fisherman switches to Rod A to enable feeding the shark.")

    # Phase 2 & 3: Power-up the shark and clear the pond.
    # We use an event-driven model.
    
    # Initialize the event queue
    # (time, event_name, description)
    events = []
    
    # Fisherman's next action is to feed the fish he just caught.
    # Feeding doesn't stop fishing, so a catch is also scheduled.
    events.append((time + 2, 'feed_shark', 'Fisherman feeds the shark.'))
    events.append((time + 7, 'catch_A', 'Fisherman catches a goldfish with Rod A.'))
    
    # Shark's natural eating schedule
    for i in range(1, 11): # Schedule natural eating every 10 mins
        events.append((i * 10, 'shark_eat_natural', 'Shark eats a free-swimming goldfish.'))

    # Main simulation loop
    while goldfish_in_pond > 0:
        if not events:
            break

        # Get the next event
        events.sort()
        current_time, event_type, description = events.pop(0)
        
        # Ignore events that are now impossible (e.g., shark eating from an empty pond)
        if 'eat' in event_type and goldfish_in_pond <= 0:
            continue

        time = current_time
        
        # Process the event
        if event_type == 'catch_A':
            goldfish_in_pond -= 1
            goldfish_in_basket += 1
            print(f"Time {time:2}: {description} {goldfish_in_pond} goldfish left in pond, {goldfish_in_basket} in basket.")
            # Immediately schedule the next actions: feeding this new fish and the next catch
            events.append((time + 2, 'feed_shark', 'Fisherman feeds the shark.'))
            events.append((time + 7, 'catch_A', 'Fisherman catches a goldfish with Rod A.'))

        elif event_type == 'feed_shark':
            if goldfish_in_basket > 0:
                goldfish_in_basket -= 1
                shark_eaten_count += 1
                shark_last_eaten_time = time
                shark_eat_times.append(time)
                print(f"Time {time:2}: {description} Shark has now eaten {shark_eaten_count} fish.")
                if shark_eaten_count >= 6 and not shark_is_strong:
                    shark_is_strong = True
                    print(f"Time {time:2}: *** Shark has eaten 6 fish and is now stronger! It will eat every 2 minutes. ***")
                    # Remove all old 'shark_eat_natural' events
                    events = [e for e in events if e[1] != 'shark_eat_natural']
                    # Schedule the next strong eat
                    events.append((time + 2, 'shark_eat_strong', 'Strong shark eats a free-swimming goldfish.'))

        elif event_type == 'shark_eat_natural':
            if not shark_is_strong:
                goldfish_in_pond -= 1
                shark_eaten_count += 1
                shark_last_eaten_time = time
                shark_eat_times.append(time)
                print(f"Time {time:2}: {description} {goldfish_in_pond} goldfish left in pond. Shark has eaten {shark_eaten_count} fish.")
                if shark_eaten_count >= 6 and not shark_is_strong:
                    shark_is_strong = True
                    print(f"Time {time:2}: *** Shark has eaten 6 fish and is now stronger! It will eat every 2 minutes. ***")
                    # Remove all old 'shark_eat_natural' events
                    events = [e for e in events if e[1] != 'shark_eat_natural']
                    # Schedule the next strong eat
                    events.append((time + 2, 'shark_eat_strong', 'Strong shark eats a free-swimming goldfish.'))

        elif event_type == 'shark_eat_strong':
            goldfish_in_pond -= 1
            shark_eaten_count += 1
            shark_last_eaten_time = time
            shark_eat_times.append(time)
            print(f"Time {time:2}: {description} {goldfish_in_pond} goldfish left in pond. Shark has eaten {shark_eaten_count} fish.")
            if goldfish_in_pond > 0:
                events.append((time + 2, 'shark_eat_strong', 'Strong shark eats a free-swimming goldfish.'))
    
    goldfish_gone_time = time
    print("-" * 60)
    print(f"All goldfish are gone from the pond at Time {goldfish_gone_time}.")
    
    # Phase 4: Calculate shark starvation
    print("\nNow, we calculate when the shark starves.")
    
    eaten_last_60_min = len([t for t in shark_eat_times if goldfish_gone_time - t <= 60])
    
    print(f"The shark's last meal was at Time {shark_last_eaten_time}.")
    print(f"In the 60 minutes before Time {goldfish_gone_time}, the shark ate {eaten_last_60_min} fish.")

    starvation_duration = 11
    extension = 4
    if eaten_last_60_min > 4:
        print(f"Since {eaten_last_60_min} is greater than 4, its survival time is extended by {extension} minutes.")
        starvation_duration += extension
        final_time = shark_last_eaten_time + starvation_duration
        print(f"The standard starvation time is 11 minutes. The total starvation time is 11 + {extension} = {starvation_duration} minutes.")
    else:
        final_time = shark_last_eaten_time + starvation_duration
        print(f"Since {eaten_last_60_min} is not greater than 4, the standard starvation time of {starvation_duration} minutes applies.")

    print("\n--- Final Calculation ---")
    print(f"The shark starves at: Last Meal Time + Starvation Duration")
    if eaten_last_60_min > 4:
        print(f"Final Time = {shark_last_eaten_time} + (11 + {extension}) = {final_time} minutes.")
    else:
        print(f"Final Time = {shark_last_eaten_time} + {starvation_duration} = {final_time} minutes.")
    
    print("\nTherefore, the earliest possible time when there are no more fish left in the pond is:")
    print(f"<<<{final_time}>>>")

if __name__ == '__main__':
    solve_fish_puzzle()