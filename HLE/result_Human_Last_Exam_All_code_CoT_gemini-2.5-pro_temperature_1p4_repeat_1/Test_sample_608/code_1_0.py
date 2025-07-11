import random

def run_monty_hall_simulations(num_simulations):
    """
    Simulates the Monty Hall problem for two strategies: staying and switching.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # Doors are represented by numbers 0, 1, 2
        doors = [0, 1, 2]
        
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # The player makes their initial choice
        player_choice = random.choice(doors)
        
        # The host opens a door that is not the player's choice and not the prize door
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)
        
        # The player decides to switch to the other unopened door
        # The switched door is the one that is not the initial choice and not the one the host opened
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # Tally wins for the "stay" strategy
        if player_choice == prize_door:
            stay_wins += 1
            
        # Tally wins for the "switch" strategy
        if switch_choice == prize_door:
            switch_wins += 1

    # --- Print the explanation and results ---
    print(f"Running {num_simulations} game simulations...\n")

    print("--- Strategy 1: Always Stay ---")
    print("You win by staying only if your initial guess was correct.")
    print("The probability of this is 1 out of 3 doors.")
    print("Theoretical Probability Equation: P(Win|Stay) = 1 / 3")
    print(f"Simulation Result: Won {stay_wins} out of {num_simulations} games.")
    print(f"Empirical Win Rate: {stay_wins / num_simulations:.2%}\n")

    print("--- Strategy 2: Always Switch ---")
    print("You win by switching if your initial guess was wrong.")
    print("The probability of this is 2 out of 3 doors.")
    print("Theoretical Probability Equation: P(Win|Switch) = 2 / 3")
    print(f"Simulation Result: Won {switch_wins} out of {num_simulations} games.")
    print(f"Empirical Win Rate: {switch_wins / num_simulations:.2%}\n")

    print("Conclusion: Yes, switching your choice roughly doubles the probability of winning.")


if __name__ == '__main__':
    run_monty_hall_simulations(10000)