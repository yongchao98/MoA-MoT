import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of trials.
    Returns the number of wins for staying vs. switching.
    """
    stay_wins = 0
    switch_wins = 0
    doors = [0, 1, 2]  # Representing the three doors

    for _ in range(num_simulations):
        # 1. Randomly place the prize behind one door
        prize_door = random.choice(doors)

        # 2. The contestant makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door.
        # The host must open a door that is NOT the player's choice and NOT the prize door.
        doors_host_can_open = [d for d in doors if d != player_choice and d != prize_door]
        opened_door = random.choice(doors_host_can_open)

        # 4. The contestant decides whether to stay or switch.
        # The switched door is the one that is not the player's choice and not the opened door.
        switched_choice = [d for d in doors if d != player_choice and d != opened_door][0]

        # Count the wins for the "stay" strategy
        if player_choice == prize_door:
            stay_wins += 1

        # Count the wins for the "switch" strategy
        if switched_choice == prize_door:
            switch_wins += 1

    return stay_wins, switch_wins

# --- Main Program ---
if __name__ == "__main__":
    simulations_to_run = 100000
    stay_wins_count, switch_wins_count = run_monty_hall_simulation(simulations_to_run)

    print(f"Running {simulations_to_run} simulations of the Monty Hall problem...")
    print("-" * 50)
    print(f"Wins if you ALWAYS STAY with your first choice: {stay_wins_count}")
    print(f"Wins if you ALWAYS SWITCH to the other door: {switch_wins_count}\n")

    print("The theoretical probabilities are:")
    print("Probability of winning by staying = 1 / 3")
    print("Probability of winning by switching = 2 / 3")
    print("\nConclusion: Yes, switching doors increases your probability of winning.")