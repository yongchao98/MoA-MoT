import random

def run_monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0
    
    doors = [0, 1, 2]

    for _ in range(num_trials):
        # 1. Place the prize behind a random door
        prize_door = random.choice(doors)
        
        # 2. Player makes an initial choice
        player_choice = random.choice(doors)
        
        # 3. Host opens a door that is not the player's choice and not the prize door
        # Create a list of doors the host can open
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)
        
        # 4. Determine the 'switch' choice. It's the door that is not the player's
        # original choice and not the one the host opened.
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # 5. Tally the wins for both strategies
        # 'Stay' strategy wins if the initial choice was the prize
        if player_choice == prize_door:
            stay_wins += 1
        
        # 'Switch' strategy wins if the switch choice was the prize
        if switch_choice == prize_door:
            switch_wins += 1

    # Calculate probabilities
    prob_stay = stay_wins / num_trials
    prob_switch = switch_wins / num_trials

    print(f"Running {num_trials} simulations of the Monty Hall problem...")
    print("-" * 50)
    print(f"Strategy: Always Stay")
    print(f"Wins: {stay_wins} out of {num_trials}")
    print(f"Winning probability by staying is {stay_wins} / {num_trials} = {prob_stay:.4f} (approximately 1/3)")
    print("-" * 50)
    print(f"Strategy: Always Switch")
    print(f"Wins: {switch_wins} out of {num_trials}")
    print(f"Winning probability by switching is {switch_wins} / {num_trials} = {prob_switch:.4f} (approximately 2/3)")


if __name__ == '__main__':
    # Run the simulation with 100,000 trials for a clear result
    run_monty_hall_simulation(100000)

<<<Yes, changing your choice increases the probability of winning from 1/3 to 2/3.>>>