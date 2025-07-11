import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of iterations.
    Calculates and prints the win percentages for both switching and staying.
    """
    doors = [1, 2, 3]
    switch_wins = 0
    stay_wins = 0

    print(f"Running {num_simulations} simulations...")

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)

        # 2. The contestant makes an initial choice
        contestant_choice = random.choice(doors)

        # 3. The host opens a door that is not the prize and not the contestant's choice
        doors_host_can_open = [d for d in doors if d != prize_door and d != contestant_choice]
        door_opened_by_host = random.choice(doors_host_can_open)

        # --- Evaluate the "Stay" Strategy ---
        if contestant_choice == prize_door:
            stay_wins += 1

        # --- Evaluate the "Switch" Strategy ---
        # The contestant switches to the only remaining closed door
        switched_choice_options = [d for d in doors if d != contestant_choice and d != door_opened_by_host]
        switched_choice = switched_choice_options[0]

        if switched_choice == prize_door:
            switch_wins += 1

    # --- Calculate and Print Results ---
    stay_win_percentage = (stay_wins / num_simulations) * 100
    switch_win_percentage = (switch_wins / num_simulations) * 100

    print("\n--- Results ---")
    print(f"Strategy: Stay with original choice")
    print(f"Wins: {stay_wins} / {num_simulations}")
    print(f"Final Equation: ({stay_wins} / {num_simulations}) * 100 = {stay_win_percentage:.2f}%")

    print(f"\nStrategy: Switch to the other door")
    print(f"Wins: {switch_wins} / {num_simulations}")
    print(f"Final Equation: ({switch_wins} / {num_simulations}) * 100 = {switch_win_percentage:.2f}%")

    print("\nConclusion: Switching doors approximately doubles your chances of winning.")


if __name__ == '__main__':
    # You can change this number to see how the probabilities hold up
    number_of_games = 10000
    run_monty_hall_simulation(number_of_games)

<<<Yes>>>