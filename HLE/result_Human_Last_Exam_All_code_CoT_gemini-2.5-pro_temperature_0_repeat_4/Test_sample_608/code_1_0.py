import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of trials.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # The doors are represented by numbers 0, 1, 2
        # Randomly place the car behind one door
        car_door = random.randint(0, 2)

        # The player makes an initial choice
        player_choice = random.randint(0, 2)

        # The host opens a door that is not the player's choice and not the car
        host_opened_door = -1
        for i in range(3):
            if i != player_choice and i != car_door:
                host_opened_door = i
                break

        # --- Strategy 1: Stay with the initial choice ---
        if player_choice == car_door:
            stay_wins += 1

        # --- Strategy 2: Switch to the other unopened door ---
        # The switched door is the one that is not the player's choice and not the one the host opened
        switch_choice = -1
        for i in range(3):
            if i != player_choice and i != host_opened_door:
                switch_choice = i
                break
        
        if switch_choice == car_door:
            switch_wins += 1

    # --- Print the results ---
    print(f"Running {num_simulations} simulations of the Monty Hall game...\n")

    print("Results for the 'Stay' strategy:")
    stay_probability = stay_wins / num_simulations
    print(f"Wins = {stay_wins}")
    print(f"Total simulations = {num_simulations}")
    print(f"Winning Probability Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}\n")

    print("Results for the 'Switch' strategy:")
    switch_probability = switch_wins / num_simulations
    print(f"Wins = {switch_wins}")
    print(f"Total simulations = {num_simulations}")
    print(f"Winning Probability Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}\n")

    print("Conclusion: The simulation shows that switching doors roughly doubles your chances of winning.")

# You can change the number of simulations to see how the probabilities converge
run_monty_hall_simulation(10000)
<<<Yes>>>