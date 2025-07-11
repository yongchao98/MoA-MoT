import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stick_wins = 0
    switch_wins = 0
    doors = [0, 1, 2]  # Representing the three doors

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)
        
        # 2. The contestant makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door that is not the prize and not the player's choice.
        # We can simplify the logic:
        # - If the player's initial choice was the prize door, switching will always lose.
        # - If the player's initial choice was a goat door, the host must open the other goat door.
        #   Switching will therefore always lead to the prize.

        # Check the "stick" strategy
        if player_choice == prize_door:
            stick_wins += 1
        
        # Check the "switch" strategy
        # You win by switching if your original choice was NOT the prize door.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate and print the probabilities
    stick_win_percentage = (stick_wins / num_simulations) * 100
    switch_win_percentage = (switch_wins / num_simulations) * 100

    print(f"Running {num_simulations} simulations of the Monty Hall problem...\n")
    print("Strategy 1: Sticking with the original choice.")
    print(f"Wins: {stick_wins} / {num_simulations}")
    print(f"Probability of winning: {stick_win_percentage:.2f}%\n")
    
    print("Strategy 2: Switching to the other door.")
    print(f"Wins: {switch_wins} / {num_simulations}")
    print(f"Probability of winning: {switch_win_percentage:.2f}%\n")

    print("Conclusion: Switching doors is the superior strategy.")

# Run the simulation with a large number of trials
monty_hall_simulation(100000)