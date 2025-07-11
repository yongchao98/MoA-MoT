import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0
    
    # A set of 3 doors
    all_doors = {1, 2, 3}
    
    for _ in range(num_simulations):
        # Step 1: Set up the game
        # Randomly place the prize behind one door
        prize_door = random.choice(list(all_doors))
        
        # Player makes their initial choice
        player_choice = random.choice(list(all_doors))
        
        # Step 2: The host opens a door
        # The host must open a door that is NOT the player's choice and NOT the prize door.
        doors_host_can_open = all_doors - {player_choice, prize_door}
        door_opened_by_host = random.choice(list(doors_host_can_open))
        
        # Step 3: Player decides to switch or stay
        # The "switch" choice is the remaining door that is not the initial choice and not the opened one.
        switched_choice = (all_doors - {player_choice, door_opened_by_host}).pop()
        
        # Step 4: Record the results
        # "Stay" strategy wins if the initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1
            
        # "Switch" strategy wins if the switched choice was the prize door
        if switched_choice == prize_door:
            switch_wins += 1
            
    # Calculate probabilities
    stay_prob = stay_wins / num_simulations
    switch_prob = switch_wins / num_simulations
    
    print(f"Simulating the game {num_simulations} times...\n")
    print("Strategy 1: Always Stay with the initial choice.")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Equation: {stay_wins} / {num_simulations} = {stay_prob:.4f}")
    print(f"This is approximately a 1/3 chance of winning.\n")
    
    print("Strategy 2: Always Switch to the other door.")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Equation: {switch_wins} / {num_simulations} = {switch_prob:.4f}")
    print(f"This is approximately a 2/3 chance of winning.\n")

    print("Conclusion: Switching doors roughly doubles your chances of winning the prize.")

# Run the simulation with 10,000 trials
monty_hall_simulation(10000)
<<<Yes>>>