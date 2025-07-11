import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to compare 'stay' vs 'switch' strategies.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0
    
    # The doors are numbered 0, 1, 2
    doors = [0, 1, 2]

    for _ in range(num_simulations):
        # 1. Place the prize behind a random door
        prize_door = random.choice(doors)
        
        # 2. The player makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door.
        # The host must open a door that is NOT the player's choice AND NOT the prize door.
        doors_host_can_open = []
        for door in doors:
            if door != player_choice and door != prize_door:
                doors_host_can_open.append(door)
        
        host_opened_door = random.choice(doors_host_can_open)

        # 4. The player decides to switch.
        # The switched choice is the door that is not the player's original choice and not the one the host opened.
        switched_choice = -1 # Initialize with an invalid value
        for door in doors:
            if door != player_choice and door != host_opened_door:
                switched_choice = door
                break

        # 5. Check the outcomes for both strategies.
        # Strategy 1: Stay with the original choice
        if player_choice == prize_door:
            stay_wins += 1
            
        # Strategy 2: Switch to the other door
        if switched_choice == prize_door:
            switch_wins += 1

    # Calculate and print the results
    stay_win_probability = stay_wins / num_simulations
    switch_win_probability = switch_wins / num_simulations

    print("Monty Hall Problem Simulation Results")
    print(f"Total simulations: {num_simulations}")
    print("-" * 40)
    print("Strategy 1: Always Stay")
    print(f"Final Equation: {stay_wins} wins / {num_simulations} simulations = {stay_win_probability:.2%} win rate")
    print("\nStrategy 2: Always Switch")
    print(f"Final Equation: {switch_wins} wins / {num_simulations} simulations = {switch_win_probability:.2%} win rate")
    print("-" * 40)
    print("\nConclusion: Switching doors approximately doubles your chances of winning.")

# Run the simulation with 10,000 trials
run_monty_hall_simulation(10000)
