import random

def monty_hall_simulation():
    """
    This function simulates the Monty Hall problem to demonstrate
    the probabilities of winning by staying versus switching.
    """
    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")

    for _ in range(num_simulations):
        # All doors are initially closed. Doors are numbered 0, 1, 2.
        all_doors = {0, 1, 2}

        # The prize is placed behind one random door.
        prize_door = random.choice(list(all_doors))

        # The player makes an initial choice.
        player_choice = random.choice(list(all_doors))

        # The host must open a door that is not the player's choice and not the prize door.
        # This is the set of doors available for the host to open.
        host_options = all_doors - {player_choice, prize_door}
        host_opens = random.choice(list(host_options))

        # The door the player can switch to is the one that is not their
        # initial choice and not the one the host opened.
        switch_choice = (all_doors - {player_choice, host_opens}).pop()

        # Tally the wins for each strategy.
        # 'Stay' strategy wins if the initial choice had the prize.
        if player_choice == prize_door:
            stay_wins += 1

        # 'Switch' strategy wins if the switched-to door has the prize.
        if switch_choice == prize_door:
            switch_wins += 1

    # Calculate and print the probabilities based on the simulation
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- Results ---")
    print(f"Strategy 1: Always Stay")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Probability of winning by staying: {stay_wins}/{num_simulations} = {stay_probability:.2%}")

    print("\nStrategy 2: Always Switch")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Probability of winning by switching: {switch_wins}/{num_simulations} = {switch_probability:.2%}")
    print("\nConclusion: Switching doors roughly doubles your chances of winning.")

if __name__ == '__main__':
    monty_hall_simulation()

<<<Yes, changing doors increases the probability of winning from 1/3 to 2/3.>>>