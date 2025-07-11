import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0
    
    # The doors are represented by numbers 0, 1, and 2
    doors = [0, 1, 2]

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)
        
        # 2. The player makes an initial choice
        player_choice = random.choice(doors)
        
        # 3. The host opens a door that is not the player's choice and not the prize door
        remaining_doors = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(remaining_doors)
        
        # 4. The player decides to switch or stay
        
        # Calculate win for the "stay" strategy
        if player_choice == prize_door:
            stay_wins += 1
            
        # Calculate win for the "switch" strategy
        # The player switches to the door that is not their choice and not the one the host opened
        switched_choice = [door for door in doors if door != player_choice and door != host_opens][0]
        if switched_choice == prize_door:
            switch_wins += 1
            
    # Calculate and print the probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations
    
    print(f"Simulating the Monty Hall game {num_simulations} times...\n")
    
    print("--- Strategy 1: Always Stay ---")
    print(f"Number of wins: {stay_wins}")
    print(f"Probability equation: {stay_wins} / {num_simulations}")
    print(f"Winning probability: {stay_probability:.4f} (approximately 1/3)\n")
    
    print("--- Strategy 2: Always Switch ---")
    print(f"Number of wins: {switch_wins}")
    print(f"Probability equation: {switch_wins} / {num_simulations}")
    print(f"Winning probability: {switch_probability:.4f} (approximately 2/3)\n")
    
    print("Conclusion: The simulation shows that switching doors significantly increases your probability of winning.")

# Run the simulation with 10,000 games
if __name__ == "__main__":
    run_monty_hall_simulation(10000)

<<<Yes>>>