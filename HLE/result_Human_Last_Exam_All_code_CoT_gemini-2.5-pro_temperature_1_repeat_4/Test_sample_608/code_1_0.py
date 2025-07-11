import random

def run_monty_hall_simulation():
    """
    This function simulates the Monty Hall problem to determine
    if switching doors increases the probability of winning.
    """
    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0
    
    print(f"--- Running Monty Hall Simulation for {num_simulations} trials ---")

    # The doors are numbered 1, 2, and 3
    doors = [1, 2, 3]

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)

        # 2. The contestant makes an initial choice
        contestant_choice = random.choice(doors)

        # 3. The host opens a door that is not the prize and not the contestant's choice
        doors_host_can_open = [door for door in doors if door != contestant_choice and door != prize_door]
        host_opens = random.choice(doors_host_can_open)

        # 4. The contestant decides to switch to the other unopened door
        switched_choice = [door for door in doors if door != contestant_choice and door != host_opens][0]

        # 5. Record the outcome for both strategies
        # 'Stay' strategy wins if the initial choice was correct
        if contestant_choice == prize_door:
            stay_wins += 1
        
        # 'Switch' strategy wins if the new choice is correct
        if switched_choice == prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    # Print the results
    print("\n--- Results ---")
    print("Strategy 1: Always Stay with the initial choice")
    print(f"Wins: {stay_wins} out of {num_simulations} trials.")
    print(f"Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}")
    print(f"Win Probability: {stay_probability:.2%}\n")

    print("Strategy 2: Always Switch to the other door")
    print(f"Wins: {switch_wins} out of {num_simulations} trials.")
    print(f"Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    print(f"Win Probability: {switch_probability:.2%}\n")
    
    print("Conclusion:")
    print("The simulation shows that switching doors approximately doubles your chances of winning the prize.")


if __name__ == '__main__':
    run_monty_hall_simulation()

<<<Yes>>>