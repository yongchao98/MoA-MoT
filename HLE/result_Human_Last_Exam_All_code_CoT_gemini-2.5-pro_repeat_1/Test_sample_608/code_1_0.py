import random

def monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The doors are numbered 1, 2, 3
        doors = [1, 2, 3]
        
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # Player makes an initial choice
        player_choice = random.choice(doors)
        
        # Host opens a door that is not the player's choice and not the prize door
        remaining_doors = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(remaining_doors)
        
        # The door the player can switch to is the one that is not their initial choice
        # and not the one the host opened.
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]
        
        # --- Check the outcome for both strategies ---
        
        # Strategy 1: Stay with the initial choice
        if player_choice == prize_door:
            stay_wins += 1
            
        # Strategy 2: Switch to the other unopened door
        if switch_choice == prize_door:
            switch_wins += 1

    # Calculate probabilities
    prob_stay = stay_wins / num_simulations
    prob_switch = switch_wins / num_simulations

    print(f"Simulating the game {num_simulations} times...\n")
    print("--- Results ---")
    print(f"Wins if you STAY with your initial choice: {stay_wins}")
    print(f"Wins if you SWITCH to the other door: {switch_wins}\n")
    
    print("--- Probabilities ---")
    print(f"Probability of winning if you STAY: {stay_wins} / {num_simulations} = {prob_stay:.4f} (or approximately 1/3)")
    print(f"Probability of winning if you SWITCH: {switch_wins} / {num_simulations} = {prob_switch:.4f} (or approximately 2/3)\n")
    
    print("Conclusion: Yes, switching doors significantly increases your probability of winning from about 33.3% to about 66.7%.")
    print("The initial choice has a 1/3 chance of being correct. The other two doors combined have a 2/3 chance. When the host reveals a goat, that 2/3 probability is concentrated on the single remaining door you can switch to.")

# Run the simulation with 10,000 trials
monty_hall_simulation(10000)
<<<Yes>>>