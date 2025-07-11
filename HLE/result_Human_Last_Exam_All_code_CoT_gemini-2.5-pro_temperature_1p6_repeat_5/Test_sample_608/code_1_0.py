import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to determine the best strategy.
    """
    stay_wins = 0
    switch_wins = 0
    
    # Run the simulation for the specified number of times
    for _ in range(num_simulations):
        # The doors are numbered 0, 1, 2
        doors = [0, 1, 2]
        
        # Randomly place the car behind one door
        car_door = random.choice(doors)
        
        # Player makes an initial random choice
        player_choice = random.choice(doors)
        
        # Host opens a door that is not the player's choice and not the car door
        # Create a list of doors the host can open
        host_options = [door for door in doors if door != player_choice and door != car_door]
        host_opens = random.choice(host_options)
        
        # Player is offered to switch to the remaining unopened door
        # The remaining door is the one that is not the player's choice and not the one the host opened
        switched_choice = [door for door in doors if door != player_choice and door != host_opens][0]
        
        # --- Tally the results for both strategies ---
        
        # Strategy 1: Player stays with their original choice
        if player_choice == car_door:
            stay_wins += 1
            
        # Strategy 2: Player switches to the other door
        if switched_choice == car_door:
            switch_wins += 1

    # --- Print the results ---
    print(f"Running Monty Hall simulation for {num_simulations} games...")
    print("-" * 50)
    
    # Results for the "Stay" strategy
    stay_probability = stay_wins / num_simulations
    print("Strategy 1: ALWAYS STAY with the initial choice.")
    print(f"Wins: {stay_wins}")
    print(f"Win Probability calculation: {stay_wins} / {num_simulations} = {stay_probability:.4f} (or {stay_probability:.2%})")

    print("-" * 50)

    # Results for the "Switch" strategy
    switch_probability = switch_wins / num_simulations
    print("Strategy 2: ALWAYS SWITCH to the other door.")
    print(f"Wins: {switch_wins}")
    print(f"Win Probability calculation: {switch_wins} / {num_simulations} = {switch_probability:.4f} (or {switch_probability:.2%})")
    
    print("-" * 50)
    print("\nConclusion: The simulation shows that switching doors approximately doubles your chances of winning the prize.")

# Run the simulation with 100,000 trials for a clear result
run_monty_hall_simulation(100000)