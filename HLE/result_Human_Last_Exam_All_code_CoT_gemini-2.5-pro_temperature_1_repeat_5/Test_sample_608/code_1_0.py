import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine if switching doors is a better strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0
    
    print(f"Running {num_simulations} simulations of the Monty Hall problem...")
    print("-" * 50)

    for _ in range(num_simulations):
        # The doors are numbered 0, 1, 2
        doors = [0, 1, 2]
        
        # Step 1: Place the prize behind a random door
        prize_door = random.choice(doors)
        
        # Step 2: The player makes an initial choice
        player_choice = random.choice(doors)
        
        # Step 3: The host opens a door that is not the player's choice and not the prize door
        # Create a list of doors the host can open
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)
        
        # Step 4: Determine the door the player would switch to
        # This is the door that is not the player's choice and not the one the host opened
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # Tally the results for both strategies
        # 'Stay' strategy wins if the initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1
            
        # 'Switch' strategy wins if the switched-to door has the prize
        if switch_choice == prize_door:
            switch_wins += 1
            
    # Calculate and print the results
    print("Strategy: Always Stay with the initial choice")
    stay_probability = stay_wins / num_simulations
    print(f"Wins: {stay_wins} out of {num_simulations} simulations.")
    print(f"Probability of winning by staying: {stay_wins} / {num_simulations} = {stay_probability:.4f} or {stay_probability:.2%}")
    print("\n")
    
    print("Strategy: Always Switch to the other door")
    switch_probability = switch_wins / num_simulations
    print(f"Wins: {switch_wins} out of {num_simulations} simulations.")
    print(f"Probability of winning by switching: {switch_wins} / {num_simulations} = {switch_probability:.4f} or {switch_probability:.2%}")
    print("-" * 50)
    
    print("\nConclusion: Yes, switching doors significantly increases your probability of winning.")
    print("The probability moves from your initial 1/3 chance to the much better 2/3 chance.")


if __name__ == '__main__':
    # Run the simulation 100,000 times for a statistically significant result
    run_monty_hall_simulation(100000)

<<<Yes>>>