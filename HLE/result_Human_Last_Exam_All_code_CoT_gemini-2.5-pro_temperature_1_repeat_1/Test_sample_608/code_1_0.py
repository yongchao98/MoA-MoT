import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0

    # The doors are represented by numbers 0, 1, 2
    doors = [0, 1, 2]

    for _ in range(num_simulations):
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)

        # Player makes an initial choice
        player_choice = random.choice(doors)

        # --- Strategy 1: Stay with the initial choice ---
        # The player wins if they initially picked the prize door.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Strategy 2: Switch to the other unopened door ---
        # The host opens a door that is not the player's choice and not the prize door.
        # If the player's initial choice was the prize door, switching will always lose.
        # If the player's initial choice was a goat door, the host must open the other goat door,
        # leaving the prize door as the only option to switch to. So, switching will always win.
        # Therefore, the "switch" strategy wins if the player's initial choice was NOT the prize door.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate and print the probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print(f"Number of simulations: {num_simulations}")
    print("\n--- Results ---")
    print(f"Winning by staying with the initial choice: {stay_wins} / {num_simulations} = {stay_probability:.4f} (or {stay_probability:.2%})")
    print(f"Winning by switching to the other door: {switch_wins} / {num_simulations} = {switch_probability:.4f} (or {switch_probability:.2%})")
    print("\nConclusion: The simulation shows that switching doors approximately doubles your chances of winning.")

# Run the simulation with a large number of trials
if __name__ == "__main__":
    monty_hall_simulation(100000)