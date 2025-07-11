import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0
    
    # Doors are represented by numbers 1, 2, and 3
    doors = [1, 2, 3]
    
    for _ in range(num_simulations):
        # 1. Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # 2. Player makes an initial choice
        player_choice = random.choice(doors)
        
        # 3. Host opens a door that is not the prize door and not the player's choice
        remaining_doors = [door for door in doors if door != player_choice and door != prize_door]
        door_opened_by_host = random.choice(remaining_doors)
        
        # 4. Determine the door the player would switch to
        # This is the door that is not the player's choice and not the one the host opened.
        switch_choice = [door for door in doors if door != player_choice and door != door_opened_by_host][0]
        
        # 5. Check the outcome for both strategies
        # 'Stay' strategy wins if the initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1
            
        # 'Switch' strategy wins if the switched-to door has the prize
        if switch_choice == prize_door:
            switch_wins += 1
            
    # Print the results
    print(f"Yes, changing doors increases your probability of winning.")
    print(f"To prove this, we simulated the game {num_simulations} times:\n")
    
    # Staying with the initial choice
    stay_probability = stay_wins / num_simulations
    print(f"Strategy 1: 'Staying' with the initial choice.")
    print(f"   - Wins: {stay_wins} out of {num_simulations}")
    print(f"   - Win probability: {stay_wins} / {num_simulations} = {stay_probability:.2%}")
    print("-" * 30)

    # Switching to the other door
    switch_probability = switch_wins / num_simulations
    print(f"Strategy 2: 'Switching' to the other unopened door.")
    print(f"   - Wins: {switch_wins} out of {num_simulations}")
    print(f"   - Win probability: {switch_wins} / {num_simulations} = {switch_probability:.2%}\n")

    print("Conclusion: The 'Switching' strategy wins about 2/3 of the time, while the 'Staying' strategy wins only about 1/3 of the time.")


# Run the simulation with 10,000 games
monty_hall_simulation(10000)