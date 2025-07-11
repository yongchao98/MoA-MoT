import random

def run_monty_hall_simulation(num_simulations=10000):
    """
    Simulates the Monty Hall problem to demonstrate whether switching
    doors increases the probability of winning.
    """
    
    print(f"Hello! This script will simulate the Monty Hall problem {num_simulations} times.")
    print("The question is: Does switching your chosen door increase your chance of winning the prize?")
    print("-" * 50)

    stick_strategy_wins = 0
    switch_strategy_wins = 0

    # The doors are represented by numbers 0, 1, and 2.
    doors = [0, 1, 2]

    for i in range(num_simulations):
        # The game master randomly places the prize behind one door.
        prize_door = random.choice(doors)

        # The player makes their initial choice.
        player_choice = random.choice(doors)

        # --- Evaluate the "sticking" strategy ---
        # If the player's initial choice was the prize door, they win by sticking.
        if player_choice == prize_door:
            stick_strategy_wins += 1

        # --- Evaluate the "switching" strategy ---
        # The core logic for switching is simpler than it seems.
        # If the player's initial choice was a goat (2/3 probability),
        # the host MUST open the other goat door. Switching to the
        # remaining door guarantees a win.
        # If the player's initial choice was the prize (1/3 probability),
        # switching guarantees a loss.
        # Therefore, the player wins by switching if their initial choice was NOT the prize.
        if player_choice != prize_door:
            switch_strategy_wins += 1

    # Calculate and display the results
    print("\nSimulation Results:\n")

    print("Strategy 1: Always STICK with the initial choice.")
    stick_win_probability = stick_strategy_wins / num_simulations
    print(f"Final Equation: {stick_strategy_wins} wins / {num_simulations} total games = {stick_win_probability:.4f}")
    print(f"This strategy wins about {stick_win_probability:.1%} of the time (approx. 1/3).\n")

    print("Strategy 2: Always SWITCH to the other door.")
    switch_win_probability = switch_strategy_wins / num_simulations
    print(f"Final Equation: {switch_strategy_wins} wins / {num_simulations} total games = {switch_win_probability:.4f}")
    print(f"This strategy wins about {switch_win_probability:.1%} of the time (approx. 2/3).\n")
    
    print("Conclusion: Yes, switching doors significantly increases your probability of winning.")


# Run the simulation when the script is executed
run_monty_hall_simulation()

<<<Yes>>>