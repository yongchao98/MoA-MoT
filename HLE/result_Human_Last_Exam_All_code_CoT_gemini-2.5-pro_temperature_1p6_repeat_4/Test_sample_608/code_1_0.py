import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # Doors are represented by numbers 1, 2, 3
        doors = [1, 2, 3]

        # Randomly place the prize behind one door
        prize_door = random.choice(doors)

        # Player makes an initial choice
        player_choice = random.choice(doors)

        # Host opens a door that is not the player's choice and not the prize door
        remaining_doors = [d for d in doors if d != player_choice and d != prize_door]
        door_opened_by_host = random.choice(remaining_doors)

        # Player is offered to switch to the other unopened door
        # The other unopened door is the one that is not the player's choice and not the one the host opened
        switched_choice = [d for d in doors if d != player_choice and d != door_opened_by_host][0]

        # --- Tally the results ---

        # "Stay" strategy: Player wins if their initial choice was the prize
        if player_choice == prize_door:
            stay_wins += 1

        # "Switch" strategy: Player wins if the switched door has the prize
        if switched_choice == prize_door:
            switch_wins += 1

    # Calculate and print the probabilities
    print(f"Number of simulations: {num_simulations}")
    print("-" * 30)

    # Staying Strategy
    stay_probability = stay_wins / num_simulations
    print("Strategy: Stay with the initial choice")
    print(f"Wins: {stay_wins}")
    print(f"Winning Probability: {stay_wins} / {num_simulations} = {stay_probability:.4f} (or {stay_probability:.2%})")
    print()

    # Switching Strategy
    switch_probability = switch_wins / num_simulations
    print("Strategy: Switch to the other door")
    print(f"Wins: {switch_wins}")
    print(f"Winning Probability: {switch_wins} / {num_simulations} = {switch_probability:.4f} (or {switch_probability:.2%})")
    print("-" * 30)
    print("\nConclusion: Switching doors approximately doubles your chance of winning.")


# Run the simulation with 10,000 games
monty_hall_simulation(10000)
