import random

def monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0
    
    # The doors are numbered 1, 2, 3
    doors = [1, 2, 3]

    print(f"Running Monty Hall simulation for {num_trials} trials...")
    print("-" * 50)

    for _ in range(num_trials):
        # 1. Game Setup: Place the prize behind a random door
        prize_door = random.choice(doors)
        
        # 2. Player's Initial Choice: Player picks a random door
        player_choice = random.choice(doors)
        
        # 3. Host's Action: Host reveals a goat door
        # The host must open a door that is NOT the player's choice and NOT the prize door.
        doors_host_can_open = []
        for d in doors:
            if d != player_choice and d != prize_door:
                doors_host_can_open.append(d)
        
        revealed_door = random.choice(doors_host_can_open)

        # 4. Determine wins for each strategy

        # "Stay" strategy: The player wins if their initial choice was the prize.
        if player_choice == prize_door:
            stay_wins += 1

        # "Switch" strategy: The player switches to the other unopened door.
        # This door is the one that is not the initial choice and not the revealed door.
        # The player wins if this switched door has the prize.
        # This is equivalent to winning if the initial choice was NOT the prize.
        if player_choice != prize_door:
            switch_wins += 1
            
    # Calculate and print the results
    stay_win_percentage = (stay_wins / num_trials) * 100
    print(f"Results for 'Always Stay' strategy:")
    print(f"Wins: {stay_wins} out of {num_trials} trials")
    print(f"Win Probability: {stay_win_percentage:.2f}% (approximately 1/3)")
    
    print("-" * 50)

    switch_win_percentage = (switch_wins / num_trials) * 100
    print(f"Results for 'Always Switch' strategy:")
    print(f"Wins: {switch_wins} out of {num_trials} trials")
    print(f"Win Probability: {switch_win_percentage:.2f}% (approximately 2/3)")

    print("-" * 50)
    print("Conclusion: Switching doors roughly doubles your chances of winning.")
    
    # Final equation as requested
    print("\nThe final probability comparison is:")
    print("2 / 3 > 1 / 3")


# Set the number of games to simulate
num_simulations = 10000
monty_hall_simulation(num_simulations)