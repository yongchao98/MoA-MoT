import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0
    
    doors = [1, 2, 3]
    
    for _ in range(num_simulations):
        # 1. Place the prize behind a random door
        prize_door = random.choice(doors)
        
        # 2. Player makes an initial choice
        player_choice = random.choice(doors)
        
        # 3. Host opens a goat door from the remaining doors
        # The host cannot open the prize door or the player's chosen door.
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)
        
        # 4. Player decides to switch or stay
        # Determine the door the player would switch to.
        # It's the one that isn't the player's choice and isn't the one the host opened.
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # Tally the wins for each strategy
        # 'Stay' strategy wins if the initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1
            
        # 'Switch' strategy wins if the switched door has the prize
        if switch_choice == prize_door:
            switch_wins += 1
            
    # Calculate and print the results
    stay_probability = (stay_wins / num_simulations) * 100
    switch_probability = (switch_wins / num_simulations) * 100
    
    print(f"Running {num_simulations} simulations of the Monty Hall Problem...\n")
    
    print("--- Strategy 1: Always Stay ---")
    print(f"Won {stay_wins} out of {num_simulations} times.")
    print(f"Win Percentage: {stay_probability:.2f}%\n")

    print("--- Strategy 2: Always Switch ---")
    print(f"Won {switch_wins} out of {num_simulations} times.")
    print(f"Win Percentage: {switch_probability:.2f}%\n")
    
    print("Conclusion: The simulation shows that switching doors approximately doubles your chances of winning.")

# Run the simulation with 10,000 trials
run_monty_hall_simulation(10000)
