import random

def run_monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_trials} simulations of the Monty Hall problem...\n")

    for _ in range(num_trials):
        # Setup the game: 0 = goat, 1 = prize
        # Let's use 0, 1, 2 for doors. Prize is behind a random door.
        prize_door = random.randint(0, 2)
        
        # Contestant makes an initial choice
        player_choice = random.randint(0, 2)
        
        # Host opens a door that is NOT the player's choice and NOT the prize door.
        host_options = [i for i in range(3) if i != player_choice and i != prize_door]
        door_opened_by_host = random.choice(host_options)
        
        # The door available to switch to is the one that isn't the player's choice
        # and isn't the one the host opened.
        switch_choice = [i for i in range(3) if i != player_choice and i != door_opened_by_host][0]
        
        # --- Evaluate the two strategies ---
        
        # Strategy 1: Stay with the initial choice
        if player_choice == prize_door:
            stay_wins += 1
            
        # Strategy 2: Switch to the other unopened door
        if switch_choice == prize_door:
            switch_wins += 1

    # Calculate and print results
    stay_probability = stay_wins / num_trials
    switch_probability = switch_wins / num_trials
    
    print("--- Results ---")
    print(f"Wins by staying with the initial choice: {stay_wins}")
    print(f"Wins by switching to the other door:   {switch_wins}")
    print("-" * 20)
    print(f"Probability of winning by staying:  {stay_probability:.4f} (approximately 1/3)")
    print(f"Probability of winning by switching: {switch_probability:.4f} (approximately 2/3)")
    print("\nConclusion: Switching doors significantly increases your chances of winning.")
    
    # Print the final equation with its numbers
    print("\nThe final equation is: P(win|switch) > P(win|stay)")
    print("Which the simulation shows is:")
    print(f"{switch_wins}/{num_trials} > {stay_wins}/{num_trials}")
    print("Or theoretically:")
    print("2/3 > 1/3")

# Run the simulation with a large number of trials for accuracy
run_monty_hall_simulation(100000)

<<<Yes>>>