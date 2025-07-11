import random

def run_monty_hall_simulation():
    """
    Simulates the Monty Hall problem to determine the best strategy
    by running the game many times and tracking the outcomes.
    """
    print("Let's analyze the Monty Hall problem by running a simulation.")
    print("We will simulate the game 10,000 times for two strategies: 'staying' and 'switching'.\n")

    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The car can be behind door 0, 1, or 2.
        car_location = random.randint(0, 2)

        # The player makes an initial choice.
        player_choice = random.randint(0, 2)

        # --- Evaluate the 'Stay' strategy ---
        # If the player's initial choice was the car, they win by staying.
        if player_choice == car_location:
            stay_wins += 1

        # --- Evaluate the 'Switch' strategy ---
        # The host opens a door that is not the player's choice and not the car.
        # If the player initially chose a goat (which happens 2/3 of the time),
        # the host must reveal the other goat. The remaining unopened door must be the car.
        # So, switching will win if the player's initial choice was NOT the car.
        if player_choice != car_location:
            switch_wins += 1

    # --- Print the results ---
    stay_probability = (stay_wins / num_simulations) * 100
    switch_probability = (switch_wins / num_simulations) * 100

    print("--- Simulation Results ---")
    print(f"Total simulations run: {num_simulations}\n")

    print("Strategy 1: Always Stay with the Initial Choice")
    print(f"Wins: {stay_wins} out of {num_simulations}")
    print(f"Winning probability equation: {stay_wins} / {num_simulations} = {stay_probability:.2f}%")
    print("\n")
    print("Strategy 2: Always Switch to the Other Door")
    print(f"Wins: {switch_wins} out of {num_simulations}")
    print(f"Winning probability equation: {switch_wins} / {num_simulations} = {switch_probability:.2f}%")
    print("\n")
    print("Conclusion: The simulation shows that switching doors roughly doubles your probability of winning.")

run_monty_hall_simulation()