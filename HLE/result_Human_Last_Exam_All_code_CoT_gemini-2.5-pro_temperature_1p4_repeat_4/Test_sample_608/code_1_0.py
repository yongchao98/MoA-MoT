import random

def run_monty_hall_simulation(num_simulations=10000):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The doors are represented by numbers 0, 1, 2
        doors = [0, 1, 2]
        
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # Player makes an initial choice
        player_choice = random.choice(doors)
        
        # Host opens a door that is not the player's choice and not the prize door
        doors_host_can_open = []
        for door in doors:
            if door != player_choice and door != prize_door:
                doors_host_can_open.append(door)
        
        host_opens = random.choice(doors_host_can_open)
        
        # --- Evaluate the "Stay" strategy ---
        if player_choice == prize_door:
            stay_wins += 1
            
        # --- Evaluate the "Switch" strategy ---
        # The player switches to the only remaining closed door
        switched_choice = -1 # initialize with a dummy value
        for door in doors:
            if door != player_choice and door != host_opens:
                switched_choice = door
                break
        
        if switched_choice == prize_door:
            switch_wins += 1

    # --- Print the results ---
    print(f"Running {num_simulations} simulations of the Monty Hall problem...")
    print("-" * 50)
    print(f"Total wins if you ALWAYS STAY: {stay_wins}")
    print(f"Total wins if you ALWAYS SWITCH: {switch_wins}")
    print("-" * 50)
    
    prob_stay = stay_wins / num_simulations
    prob_switch = switch_wins / num_simulations
    
    print(f"Calculated probability of winning by staying: {prob_stay:.4f} (Theoretical: 1/3)")
    print(f"Calculated probability of winning by switching: {prob_switch:.4f} (Theoretical: 2/3)")
    print("-" * 50)
    
    print("Conclusion: Yes, changing doors increases your probability of winning.")
    print("The final comparison shows that the probability of winning by switching is greater than by staying.")
    
    # Printing each number in the final equation as requested
    win_switch_numerator = 2
    win_switch_denominator = 3
    win_stay_numerator = 1
    win_stay_denominator = 3
    print(f"Final Equation: {win_switch_numerator}/{win_switch_denominator} > {win_stay_numerator}/{win_stay_denominator}")


# Run the simulation
run_monty_hall_simulation()