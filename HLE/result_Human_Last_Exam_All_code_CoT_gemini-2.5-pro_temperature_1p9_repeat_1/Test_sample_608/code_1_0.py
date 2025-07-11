import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem to demonstrate the winning probabilities
    of staying versus switching doors.
    """
    stay_wins = 0
    switch_wins = 0
    doors = [1, 2, 3]

    print(f"Running {num_simulations} simulations of the Monty Hall problem...\n")

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        winning_door = random.choice(doors)

        # 2. The contestant makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door that is not the winning door and not the player's choice
        doors_host_can_open = [d for d in doors if d != winning_door and d != player_choice]
        host_opens = random.choice(doors_host_can_open)

        # 4. The player decides to switch to the other unopened door
        # This is the door that is not the player's initial choice and not the one the host opened.
        switched_choice = [d for d in doors if d != player_choice and d != host_opens][0]

        # 5. Tally the wins for both strategies
        # 'Stay' strategy wins if the initial choice was the winning door.
        if player_choice == winning_door:
            stay_wins += 1

        # 'Switch' strategy wins if the switched choice was the winning door.
        if switched_choice == winning_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("--- Results ---")
    print(f"Strategy: Always Stay")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    # The final equation here is the probability calculation.
    print(f"Probability of winning: {stay_wins} / {num_simulations} = {stay_probability:.4f} (or approximately 1/3)")

    print("\n" + "="*20 + "\n")

    print(f"Strategy: Always Switch")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    # The final equation here is the probability calculation.
    print(f"Probability of winning: {switch_wins} / {num_simulations} = {switch_probability:.4f} (or approximately 2/3)")
    print("\n--- Conclusion ---")
    print("Switching doors significantly increases your probability of winning.")


if __name__ == "__main__":
    # You can change the number of simulations to see how the results converge.
    number_of_games = 10000
    run_monty_hall_simulation(number_of_games)
