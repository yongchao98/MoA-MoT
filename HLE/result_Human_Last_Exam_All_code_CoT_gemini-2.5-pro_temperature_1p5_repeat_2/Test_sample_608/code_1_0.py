import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_simulations} simulations...")

    for _ in range(num_simulations):
        # The doors are represented by numbers 0, 1, and 2
        # Randomly place the prize behind one door
        prize_door = random.randint(0, 2)
        
        # Player makes an initial choice
        player_choice = random.randint(0, 2)

        # --- Evaluate the "Stay" strategy ---
        # If the player's choice is the prize door, they win by staying.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Evaluate the "Switch" strategy ---
        # The host opens a door that is not the player's choice and not the prize.
        # If the player initially chose the prize (1/3 of the time), switching loses.
        # If the player initially chose a goat (2/3 of the time), the host must reveal
        # the other goat, so switching to the last remaining door wins.
        # Therefore, we win by switching if our initial choice was NOT the prize door.
        if player_choice != prize_door:
            switch_wins += 1
            
    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Results ---")
    print("\nStrategy 1: Always Stay")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}")
    print(f"Win Probability: {stay_probability:.2%}")

    print("\nStrategy 2: Always Switch")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    print(f"Win Probability: {switch_probability:.2%}")

    if switch_probability > stay_probability:
        print("\nConclusion: Switching doors significantly increases your probability of winning.")
    else:
        print("\nConclusion: The probabilities are about the same.")


# Run the simulation with 100,000 trials
monty_hall_simulation(100000)
