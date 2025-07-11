import random

def run_monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0
    
    # Doors are represented by numbers 0, 1, 2
    doors = [0, 1, 2]

    for _ in range(num_trials):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)
        
        # 2. The contestant makes an initial choice
        player_choice = random.choice(doors)
        
        # 3. The host opens a door that is not the prize door and not the player's choice
        # Create a list of doors the host can open
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)
        
        # 4. The player decides to stay or switch
        # The 'stay' choice is just the original player_choice
        if player_choice == prize_door:
            stay_wins += 1
            
        # The 'switch' choice is the remaining unopened door
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]
        if switch_choice == prize_door:
            switch_wins += 1

    print(f"Simulating {num_trials} trials of the Monty Hall problem...")
    
    # Calculate and print the results for the "stay" strategy
    stay_probability = stay_wins / num_trials
    print(f"Strategy: Don't Switch")
    print(f"Wins: {stay_wins} / {num_trials}")
    print(f"Probability of winning: {stay_probability:.4f}")
    
    print("-" * 30)

    # Calculate and print the results for the "switch" strategy
    switch_probability = switch_wins / num_trials
    print(f"Strategy: Switch")
    print(f"Wins: {switch_wins} / {num_trials}")
    print(f"Probability of winning: {switch_probability:.4f}")


if __name__ == '__main__':
    # Run the simulation with 100,000 trials for a clear result
    run_monty_hall_simulation(100000)
