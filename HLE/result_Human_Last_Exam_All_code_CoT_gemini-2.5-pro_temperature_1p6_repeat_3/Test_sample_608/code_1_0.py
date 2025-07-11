import random

def simulate_monty_hall(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.

    Args:
        num_simulations (int): The number of times to run the simulation.
    """
    if num_simulations <= 0:
        print("Number of simulations must be a positive integer.")
        return

    stay_wins = 0
    switch_wins = 0
    
    doors = [0, 1, 2]  # Represent the three doors

    for _ in range(num_simulations):
        # 1. Place the prize behind a random door
        prize_door = random.choice(doors)

        # 2. The player makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door.
        # The host must open a door that is NOT the player's choice and NOT the prize door.
        # This step is implicitly handled by the logic below.

        # 4. Evaluate the two strategies
        
        # Strategy 1: The player stays with their original choice.
        # They win if their initial choice was the prize door.
        if player_choice == prize_door:
            stay_wins += 1

        # Strategy 2: The player switches.
        # They win if their initial choice was NOT the prize door.
        # Because if they initially picked a goat, the host reveals the other goat,
        # leaving the prize behind the door they can switch to.
        if player_choice != prize_door:
            switch_wins += 1

    # Calculate and print the probabilities
    stay_prob = stay_wins / num_simulations
    switch_prob = switch_wins / num_simulations
    
    print(f"Running {num_simulations} simulations of the Monty Hall problem...\n")
    print("Strategy 1: Always Stay with the initial choice")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Probability of winning: {stay_wins} / {num_simulations} = {stay_prob:.4f}\n")
    
    print("Strategy 2: Always Switch to the other door")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Probability of winning: {switch_wins} / {num_simulations} = {switch_prob:.4f}\n")
    
    print("Conclusion: Switching doors roughly doubles your probability of winning.")


# You can change the number of simulations to see the effect
simulate_monty_hall(10000)
