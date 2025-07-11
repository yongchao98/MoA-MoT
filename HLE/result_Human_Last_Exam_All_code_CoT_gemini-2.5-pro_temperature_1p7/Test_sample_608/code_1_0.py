import random

def run_simulation(num_trials, switch_choice):
    """
    Runs a simulation of the Monty Hall problem.
    - num_trials: The number of times to run the simulation.
    - switch_choice: A boolean indicating if the player switches their choice.
    Returns the number of wins.
    """
    wins = 0
    # There are 3 doors, represented by numbers 0, 1, and 2.
    doors = [0, 1, 2]

    for _ in range(num_trials):
        # Step 1: The prize is placed behind a random door.
        car_door = random.choice(doors)
        
        # Step 2: The player makes an initial choice.
        player_choice = random.choice(doors)
        
        # Step 3: The host opens a door that is NOT the player's choice and NOT the car door.
        doors_host_can_open = [door for door in doors if door != player_choice and door != car_door]
        host_opens = random.choice(doors_host_can_open)

        # Step 4: The player decides to switch or stay.
        if switch_choice:
            # The player switches to the only other door that is not their original choice and not the one the host opened.
            final_choice = [door for door in doors if door != player_choice and door != host_opens][0]
        else:
            # The player stays with their original choice.
            final_choice = player_choice
            
        # Check if the final choice is the car door.
        if final_choice == car_door:
            wins += 1
            
    return wins

# Set the total number of games to simulate for each strategy.
num_simulations = 10000

print("Yes, changing your choice increases your probability of winning.\n")
print("To demonstrate this, we will simulate the game 10,000 times for each strategy.\n")

# --- STRATEGY 1: ALWAYS STAY ---
wins_staying = run_simulation(num_simulations, switch_choice=False)
prob_staying = wins_staying / num_simulations

print("--- Strategy 1: Always STAY with the initial choice ---")
print(f"Out of {num_simulations} games, staying resulted in {wins_staying} wins.")
print("The winning probability is calculated as:")
print(f"{wins_staying} / {num_simulations} = {prob_staying:.4f}")
print("This is approximately a 1/3 chance of winning.\n")

# --- STRATEGY 2: ALWAYS SWITCH ---
wins_switching = run_simulation(num_simulations, switch_choice=True)
prob_switching = wins_switching / num_simulations

print("--- Strategy 2: Always SWITCH to the other door ---")
print(f"Out of {num_simulations} games, switching resulted in {wins_switching} wins.")
print("The winning probability is calculated as:")
print(f"{wins_switching} / {num_simulations} = {prob_switching:.4f}")
print("This is approximately a 2/3 chance of winning.\n")

print("Conclusion: The simulation shows that switching doors nearly doubles your chance of winning.")
