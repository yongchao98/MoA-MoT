import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")

    for _ in range(num_simulations):
        # Doors are 0, 1, 2. The prize is randomly behind one of them.
        prize_door = random.randint(0, 2)
        
        # Contestant makes an initial choice.
        player_choice = random.randint(0, 2)

        # --- Strategy 1: Stay with the initial choice ---
        # You win if your initial choice was the prize door.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Strategy 2: Switch to the other door ---
        # The host opens a door that is not the player's choice and not the prize door.
        # The player then switches.
        # This means the player wins if their initial choice was NOT the prize door.
        # If the player initially chose a goat, the host opens the other goat door,
        # leaving the prize door as the only option to switch to.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Results ---")
    print(f"Strategy 1: Always Stay")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Winning Probability equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}")

    print(f"\nStrategy 2: Always Switch")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Winning Probability equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    
    print("\nConclusion: Yes, switching doors approximately doubles your probability of winning.")


# Run the simulation with 100,000 trials for a clear result
monty_hall_simulation(100000)

<<<Yes>>>