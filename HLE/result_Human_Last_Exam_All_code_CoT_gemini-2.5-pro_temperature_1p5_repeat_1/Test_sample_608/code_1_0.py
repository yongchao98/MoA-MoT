import random

def monty_hall_simulation():
    """
    Simulates the Monty Hall problem to determine if switching doors
    improves the probability of winning.
    """
    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0
    doors = [1, 2, 3]

    print(f"Running {num_simulations} simulations of the Monty Hall problem...")

    for _ in range(num_simulations):
        # The prize is randomly placed behind one door
        prize_door = random.choice(doors)

        # The contestant makes a random initial choice
        player_choice = random.choice(doors)

        # The host must open a door that is not the player's choice and not the prize door
        host_options = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(host_options)

        # The door available to switch to is the one that is not the player's choice and not the one the host opened
        switch_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # Calculate wins for the "Stay" strategy
        if player_choice == prize_door:
            stay_wins += 1

        # Calculate wins for the "Switch" strategy
        if switch_choice == prize_door:
            switch_wins += 1
            
    # Calculate and print the results
    stay_probability = (stay_wins / num_simulations) * 100
    switch_probability = (switch_wins / num_simulations) * 100

    print("\n--- Results ---")
    print("\nStrategy 1: Always Stay")
    print(f"Wins: {stay_wins} / {num_simulations}")
    print(f"Probability of winning by staying: {stay_probability:.2f}%")

    print("\nStrategy 2: Always Switch")
    print(f"Wins: {switch_wins} / {num_simulations}")
    print(f"Probability of winning by switching: {switch_probability:.2f}%")

    print("\nConclusion: The simulation shows that switching doors significantly increases the probability of winning.")


if __name__ == '__main__':
    monty_hall_simulation()
<<<Yes>>>