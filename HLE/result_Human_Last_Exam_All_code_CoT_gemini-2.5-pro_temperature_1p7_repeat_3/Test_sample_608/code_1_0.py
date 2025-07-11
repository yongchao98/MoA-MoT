import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of trials.

    Args:
        num_simulations (int): The number of times to run the simulation.

    Returns:
        None. Prints the results of the simulation.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The prize is behind one of three doors (1, 2, or 3)
        prize_door = random.randint(1, 3)

        # The contestant makes an initial choice
        player_choice = random.randint(1, 3)

        # --- "Stay" strategy ---
        # If the player's initial choice was the prize door, they win by staying.
        if player_choice == prize_door:
            stay_wins += 1

        # --- "Switch" strategy ---
        # The host opens a door that is not the player's choice and not the prize.
        # If the player initially chose a goat (player_choice != prize_door),
        # the host is forced to open the other goat door.
        # Switching to the last remaining door guarantees a win.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_win_probability = stay_wins / num_simulations
    switch_win_probability = switch_wins / num_simulations

    print(f"Number of simulations: {num_simulations}\n")
    print("--- Strategy: Stay with the initial choice ---")
    print(f"Wins: {stay_wins}")
    print(f"Winning Probability: {stay_win_probability:.4f} (approximately 1/3)")
    print("\n--- Strategy: Switch to the other door ---")
    print(f"Wins: {switch_wins}")
    print(f"Winning Probability: {switch_win_probability:.4f} (approximately 2/3)")

    # The final equation is the comparison of the two probabilities.
    # We will show it as a simplified fraction.
    print("\n--- Final Answer ---")
    print("Does changing your choice increase your probability of winning?")
    print("Yes. The probability of winning increases from 1/3 to 2/3.")


if __name__ == '__main__':
    # Run the simulation with 10,000 trials for a clear result
    monty_hall_simulation(10000)
