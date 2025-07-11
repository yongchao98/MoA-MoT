import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of iterations
    to demonstrate the probabilities of winning for staying vs. switching.
    """
    stay_wins = 0
    switch_wins = 0
    doors = [1, 2, 3]

    print(f"Running {num_simulations:,} simulations of the Monty Hall problem...\n")

    for _ in range(num_simulations):
        # Step 1: The prize is placed behind a random door.
        prize_door = random.choice(doors)

        # Step 2: The contestant makes an initial choice.
        player_choice = random.choice(doors)

        # --- Analyze the "Stay" Strategy ---
        # The player wins if their initial choice was the prize door.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Analyze the "Switch" Strategy ---
        # The player wins by switching if their initial choice was NOT the prize door.
        # This is because the host will always reveal the other non-prize door,
        # leaving the prize door as the only option to switch to.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    # Print the final results of the simulation
    print("--- Results ---")
    print("Strategy 1: Always Stay with the initial choice")
    print(f"Final Equation: {stay_wins} wins / {num_simulations} games = {stay_probability:.4f}")
    print(f"Win Probability: {stay_probability:.2%}\n")

    print("Strategy 2: Always Switch to the other door")
    print(f"Final Equation: {switch_wins} wins / {num_simulations} games = {switch_probability:.4f}")
    print(f"Win Probability: {switch_probability:.2%}\n")

    print("Conclusion: The simulation confirms that switching doors approximately doubles your chances of winning.")

# Set the number of games to simulate and run the function
number_of_games = 10000
monty_hall_simulation(number_of_games)