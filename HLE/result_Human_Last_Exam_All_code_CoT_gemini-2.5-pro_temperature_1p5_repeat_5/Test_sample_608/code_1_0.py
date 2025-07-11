import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The prize is behind one of the three doors (0, 1, 2)
        prize_door = random.randint(0, 2)
        
        # Contestant makes their initial choice
        player_choice = random.randint(0, 2)

        # The host must now open a door that is:
        # 1. Not the player's choice
        # 2. Not the prize door
        # This logic is implicitly handled by analyzing the outcomes.

        # --- Analyze the "Stay" Strategy ---
        # The player wins if their initial choice was the prize door.
        if player_choice == prize_door:
            stay_wins += 1

        # --- Analyze the "Switch" Strategy ---
        # If the player's initial choice was NOT the prize, the host is forced
        # to open the only other door that has a goat. The remaining unopened
        # door MUST have the prize. Therefore, switching wins.
        if player_choice != prize_door:
            switch_wins += 1

    # --- Print Results ---
    print(f"Simulating the Monty Hall game {num_simulations} times...")
    print("="*50)

    # Calculate and print the probability for the "Stay" strategy
    stay_probability = stay_wins / num_simulations
    print("Results for the 'Stay' Strategy:")
    print(f"Total Wins: {stay_wins}")
    print(f"Probability Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}")
    print(f"You have an approximate 1/3 chance of winning if you stay.")
    
    print("-"*50)

    # Calculate and print the probability for the "Switch" strategy
    switch_probability = switch_wins / num_simulations
    print("Results for the 'Switch' Strategy:")
    print(f"Total Wins: {switch_wins}")
    print(f"Probability Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    print(f"You have an approximate 2/3 chance of winning if you switch.")
    
    print("="*50)
    print("\nConclusion: Yes, switching doors significantly increases your probability of winning from ~33.3% to ~66.7%.")


# Run the simulation with 10,000 trials
run_monty_hall_simulation(10000)