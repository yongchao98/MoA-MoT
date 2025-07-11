import random

def run_monty_hall_simulation(num_simulations=100000):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    
    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    
    stay_wins = 0
    switch_wins = 0
    
    # We will loop through the number of specified simulations
    for _ in range(num_simulations):
        # Setup the game: three doors, with numbers 0, 1, 2
        doors = [0, 1, 2]
        
        # Step 1: The prize is randomly placed behind one door.
        prize_door = random.choice(doors)
        
        # Step 2: The contestant makes a random initial choice.
        player_initial_choice = random.choice(doors)
        
        # Step 3: The host opens a door.
        # The host must open a door that is NOT the player's choice AND NOT the prize door.
        host_options = [door for door in doors if door != player_initial_choice and door != prize_door]
        door_opened_by_host = random.choice(host_options)
        
        # Step 4: The player's final choice is determined for both strategies.
        # The 'stay' choice is the initial choice.
        # The 'switch' choice is the one remaining door that is not the initial choice and not the one the host opened.
        switch_choice = [door for door in doors if door != player_initial_choice and door != door_opened_by_host][0]
        
        # Tally the results for both strategies
        if player_initial_choice == prize_door:
            stay_wins += 1
            
        if switch_choice == prize_door:
            switch_wins += 1
            
    # Calculate and display the final probabilities
    print(f"Simulating the Monty Hall game for {num_simulations} rounds...\n")

    print("--- Strategy 1: Always STAY with the initial choice ---")
    stay_win_probability = stay_wins / num_simulations
    print(f"Wins: {stay_wins}")
    print(f"Probability of winning: {stay_wins} / {num_simulations} = {stay_win_probability:.4f}")
    print(f"(This is approximately 1/3)\n")
    
    print("--- Strategy 2: Always SWITCH to the other door ---")
    switch_win_probability = switch_wins / num_simulations
    print(f"Wins: {switch_wins}")
    print(f"Probability of winning: {switch_wins} / {num_simulations} = {switch_win_probability:.4f}")
    print(f"(This is approximately 2/3)\n")
    
    print("Conclusion: The simulation shows that switching your choice roughly doubles your chances of winning the prize.")

# Run the simulation
if __name__ == "__main__":
    run_monty_hall_simulation()
