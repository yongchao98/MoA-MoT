import random

def monty_hall_simulation(num_simulations=10000):
    """
    Simulates the Monty Hall problem to demonstrate the probabilities
    of winning by either staying or switching doors.
    """
    stay_wins = 0
    switch_wins = 0

    # Run the simulation for the specified number of times
    for _ in range(num_simulations):
        # Let's represent the doors with numbers 1, 2, and 3
        doors = [1, 2, 3]

        # Step 1: The prize is randomly placed behind one door
        prize_door = random.choice(doors)

        # Step 2: The player makes an initial choice
        player_choice = random.choice(doors)

        # Step 3: The host opens a door.
        # The host must open a door that is NOT the player's choice and NOT the prize door.
        doors_host_can_open = [door for door in doors if door != player_choice and door != prize_door]
        # host_opens = random.choice(doors_host_can_open) # This is not needed for the logic, but it's what happens.

        # Step 4: Determine the outcome for both strategies.

        # The 'Stay' strategy wins if the player's initial choice was the prize.
        if player_choice == prize_door:
            stay_wins += 1

        # The 'Switch' strategy wins if the player's initial choice was NOT the prize.
        # If the player initially picked a goat, the host is forced to open the other goat door.
        # Switching will then always lead to the prize door.
        if player_choice != prize_door:
            switch_wins += 1

    print(f"Running {num_simulations} simulations of the Monty Hall game...\n")
    print("Strategy 1: Always Stay")
    print(f"Wins: {stay_wins}")
    print(f"Losses: {num_simulations - stay_wins}")

    print("\nStrategy 2: Always Switch")
    print(f"Wins: {switch_wins}")
    print(f"Losses: {num_simulations - switch_wins}")

    # Calculate and display the final probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print("\n--- RESULTS ---")
    print("Yes, changing doors increases your probability of winning.")
    print("The simulation shows the following probabilities:")
    print("\nWinning probability if you STAY with your initial choice:")
    print(f"{stay_wins} / {num_simulations} = {stay_probability:.4f} (approximately 1/3)")
    print("\nWinning probability if you SWITCH to the other door:")
    print(f"{switch_wins} / {num_simulations} = {switch_probability:.4f} (approximately 2/3)")

# Execute the simulation
monty_hall_simulation()
<<<Yes>>>