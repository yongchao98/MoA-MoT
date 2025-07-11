import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for two strategies: staying and switching.
    """
    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")

    for _ in range(num_simulations):
        # The doors are represented by numbers 1, 2, and 3
        doors = [1, 2, 3]

        # Randomly place the prize behind one door
        prize_door = random.choice(doors)

        # The player makes an initial choice
        player_choice = random.choice(doors)

        # --- Strategy 1: Stay with the initial choice ---
        # If the player's initial choice was the prize door, it's a win.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Strategy 2: Switch to the other door ---
        # In this strategy, the player switches. The player wins if their
        # initial choice was NOT the prize.
        # Why? If the player initially picked a goat, the host MUST open the
        # other goat door, leaving the prize door as the only option to switch to.
        # So, a switch wins if the initial guess was wrong.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Simulation Results ---")
    print(f"Total simulations: {num_simulations}")
    print(f"Wins by staying with the initial choice: {stay_wins}")
    print(f"Wins by switching to the other door: {switch_wins}")

    print(f"\nProbability of winning if you stay: {stay_probability:.4f} (or ~1/3)")
    print(f"Probability of winning if you switch: {switch_probability:.4f} (or ~2/3)")

    print("\n--- Final Equation ---")
    print("Probability of winning by switching is roughly 2/3.")
    print("Probability of winning by staying is roughly 1/3.")
    print("\nThe numbers in the final theoretical equation are:")
    # Printing the numbers in the equation P(switch) = 2/3 > P(stay) = 1/3
    print("2")
    print("3")
    print("1")
    print("3")


if __name__ == '__main__':
    # You can change this number to see how it affects the results
    number_of_games = 10000
    run_monty_hall_simulation(number_of_games)
