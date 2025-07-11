import random

def monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0

    print(f"Running simulation for {num_trials} trials...")

    for _ in range(num_trials):
        # The doors are represented by numbers 0, 1, 2
        doors = [0, 1, 2]
        
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # The contestant makes an initial choice
        player_choice = random.choice(doors)

        # Evaluate the 'stay' strategy
        # A win occurs if the player's initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1

        # Evaluate the 'switch' strategy
        # A win occurs if the player's initial choice was NOT the prize door.
        # Why? If you initially picked a goat (which happens 2/3 of the time),
        # the host opens the other goat door, leaving the prize as the only option to switch to.
        # So, switching from a goat door is a guaranteed win.
        if player_choice != prize_door:
            switch_wins += 1

    # --- Print Results ---
    print("\n--- Simulation Results ---")

    # Stay strategy results
    stay_losses = num_trials - stay_wins
    stay_win_prob = stay_wins / num_trials
    print("\nStrategy: Always Stay")
    print(f"Wins: {stay_wins}")
    print(f"Losses: {stay_losses}")
    print(f"Equation for Win Probability: {stay_wins} / {num_trials} = {stay_win_prob:.4f}")
    
    # Switch strategy results
    switch_losses = num_trials - switch_wins
    switch_win_prob = switch_wins / num_trials
    print("\nStrategy: Always Switch")
    print(f"Wins: {switch_wins}")
    print(f"Losses: {switch_losses}")
    print(f"Equation for Win Probability: {switch_wins} / {num_trials} = {switch_win_prob:.4f}")

    print("\n--- Conclusion ---")
    print("Yes, changing your choice approximately doubles your probability of winning from ~1/3 to ~2/3.")

# Run the simulation with 10,000 trials
monty_hall_simulation(10000)
