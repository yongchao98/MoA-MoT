import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of trials.

    This simplified simulation is based on the core logic:
    - If the player's initial choice is correct (1/3 probability), staying wins and switching loses.
    - If the player's initial choice is incorrect (2/3 probability), the host's action forces the
      remaining door to be the winner, so switching wins and staying loses.
    """
    stay_wins = 0
    switch_wins = 0
    num_doors = 3

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")

    for _ in range(num_simulations):
        # We don't need to simulate all the doors, just whether the first choice was correct.
        # There's a 1 in 3 chance the player initially picked the correct door.
        # We can simulate this with a random number choice.
        if random.randint(1, num_doors) == 1:
            # Player's initial choice was the prize.
            stay_wins += 1
        else:
            # Player's initial choice was a goat.
            # Host reveals the other goat, so switching wins.
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Results ---")
    print(f"Strategy: Stay with the initial choice")
    print(f"Wins: {stay_wins} out of {num_simulations} games")
    print(f"Probability Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f} (or approximately 1/3)")

    print(f"\nStrategy: Switch to the other door")
    print(f"Wins: {switch_wins} out of {num_simulations} games")
    print(f"Probability Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f} (or approximately 2/3)")

    print("\nConclusion: The simulation shows that switching doors significantly increases the probability of winning.")


# Run the simulation with 10,000 games
run_monty_hall_simulation(10000)
