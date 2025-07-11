import random

def run_monty_hall_simulation():
    """
    This function simulates the Monty Hall problem to determine
    the best strategy: staying or switching.
    """
    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0

    print("Running a simulation of the Monty Hall Problem...")
    print(f"Number of games to simulate: {num_simulations}\n")

    for _ in range(num_simulations):
        # The doors are numbered 0, 1, 2.
        # Let's randomly place the car behind one door.
        car_location = random.randint(0, 2)

        # The player makes their initial choice.
        player_choice = random.randint(0, 2)

        # The host must open a door that is not the player's choice
        # and does not have the car.
        host_options = []
        for i in range(3):
            if i != player_choice and i != car_location:
                host_options.append(i)
        
        # The host reveals a goat by opening one of the available doors.
        host_opens = random.choice(host_options)

        # --- Strategy 1: The player stays with their original choice. ---
        if player_choice == car_location:
            stay_wins += 1

        # --- Strategy 2: The player switches to the other unopened door. ---
        # The switched door is the one not chosen by the player and not opened by the host.
        options_for_switch = [0, 1, 2]
        options_for_switch.remove(player_choice)
        options_for_switch.remove(host_opens)
        switched_choice = options_for_switch[0]

        if switched_choice == car_location:
            switch_wins += 1

    # --- Calculate and display the final results ---
    print("---------- RESULTS ----------")
    
    # Staying Strategy Results
    stay_probability = stay_wins / num_simulations
    print("Strategy: ALWAYS STAY with the initial choice")
    print(f"Wins by staying: {stay_wins}")
    print(f"Final Equation: {stay_wins} / {num_simulations} = {stay_probability:.4f}")
    print(f"This is a win probability of {stay_probability:.2%}.\n")

    # Switching Strategy Results
    switch_probability = switch_wins / num_simulations
    print("Strategy: ALWAYS SWITCH to the other door")
    print(f"Wins by switching: {switch_wins}")
    print(f"Final Equation: {switch_wins} / {num_simulations} = {switch_probability:.4f}")
    print(f"This is a win probability of {switch_probability:.2%}.\n")
    
    print("Conclusion: Yes, changing your door roughly doubles your probability of winning.")


if __name__ == '__main__':
    run_monty_hall_simulation()