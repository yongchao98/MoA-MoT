import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0
    doors = [0, 1, 2]

    for _ in range(num_simulations):
        # The prize is randomly placed behind one of the three doors.
        prize_door = random.randint(0, 2)

        # The contestant makes their initial choice.
        player_choice = random.randint(0, 2)

        # The host must open a door that is not the player's choice and not the prize door.
        doors_to_open = [door for door in doors if door != player_choice and door != prize_door]
        host_opens = random.choice(doors_to_open)

        # The contestant decides whether to stay or switch.
        # The switched choice is the door that is not the player's choice and not the one the host opened.
        switched_choice = [door for door in doors if door != player_choice and door != host_opens][0]

        # Check the outcome for the "stay" strategy.
        if player_choice == prize_door:
            stay_wins += 1

        # Check the outcome for the "switch" strategy.
        if switched_choice == prize_door:
            switch_wins += 1

    # Calculate and print the probabilities
    stay_probability = stay_wins / num_simulations
    switch_probability = switch_wins / num_simulations

    print(f"Simulating {num_simulations} games...\n")
    
    print("--- Strategy: Always Stay with the Initial Choice ---")
    print(f"Wins: {stay_wins}")
    print(f"Final Equation: {stay_wins} wins / {num_simulations} total games")
    print(f"Win Probability: {stay_probability:.4f} (approximately 1/3)\n")
    
    print("--- Strategy: Always Switch Doors ---")
    print(f"Wins: {switch_wins}")
    print(f"Final Equation: {switch_wins} wins / {num_simulations} total games")
    print(f"Win Probability: {switch_probability:.4f} (approximately 2/3)\n")
    
    print("Conclusion: Yes, switching doors significantly increases your probability of winning.")


if __name__ == '__main__':
    # You can change this number to see how the probabilities converge.
    number_of_games = 10000
    run_monty_hall_simulation(number_of_games)
