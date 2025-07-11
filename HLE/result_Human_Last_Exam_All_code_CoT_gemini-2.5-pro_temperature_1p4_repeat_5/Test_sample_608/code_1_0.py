import random

def run_monty_hall_simulation():
    """
    This function simulates the Monty Hall problem to demonstrate the
    probabilities of winning by sticking with the initial choice versus switching.
    """
    num_simulations = 10000
    stick_wins = 0
    switch_wins = 0

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")
    print("-" * 60)

    for _ in range(num_simulations):
        # Setup: Doors are indexed 0, 1, 2. Place the prize randomly.
        prize_door = random.randint(0, 2)
        
        # Player makes an initial choice
        player_choice = random.randint(0, 2)
        
        # Host opens a door.
        # The host must open a door that is not the player's choice and not the prize door.
        doors_to_open = [door for door in range(3) if door != player_choice and door != prize_door]
        host_opens = random.choice(doors_to_open)
        
        # Determine the switched door choice.
        # It's the only door left that isn't the player's choice or the host's open door.
        switch_choice = [door for door in range(3) if door != player_choice and door != host_opens][0]

        # Tally the wins for each strategy
        # Strategy 1: Stick with the original choice
        if player_choice == prize_door:
            stick_wins += 1
        
        # Strategy 2: Switch to the other unopened door
        if switch_choice == prize_door:
            switch_wins += 1

    # Calculate and print the results
    stick_probability = stick_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Strategy 1: Always Stick with the First Choice ---")
    print(f"Number of wins: {stick_wins}")
    # Output each number in the final equation
    print(f"Winning probability equation: {stick_wins} / {num_simulations} = {stick_probability:.4f}")
    print(f"This is approximately a {stick_probability:.1%} chance of winning.")

    print("\n--- Strategy 2: Always Switch to the Other Door ---")
    print(f"Number of wins: {switch_wins}")
    # Output each number in the final equation
    print(f"Winning probability equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    print(f"This is approximately a {switch_probability:.1%} chance of winning.")

    print("-" * 60)
    print("\nConclusion: The simulation shows that switching doors roughly doubles your chances of winning.")

# Run the simulation
run_monty_hall_simulation()
<<<Yes>>>