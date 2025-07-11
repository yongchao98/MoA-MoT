import random

def run_monty_hall_simulation():
    """
    Simulates the Monty Hall problem to determine if switching doors
    increases the probability of winning.
    """
    num_trials = 10000
    doors = [1, 2, 3]

    stay_wins = 0
    switch_wins = 0

    print(f"Running {num_trials} simulations of the Monty Hall problem...\n")

    for _ in range(num_trials):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)

        # 2. The contestant makes an initial choice
        contestant_choice = random.choice(doors)

        # 3. The host must open a door that is not the prize door
        #    and not the contestant's choice.
        host_options = [door for door in doors if door != prize_door and door != contestant_choice]
        host_opens = random.choice(host_options)

        # 4. Let's see the results for both staying and switching.

        # -- Result if the contestant 'stays' with the initial choice --
        if contestant_choice == prize_door:
            stay_wins += 1

        # -- Result if the contestant 'switches' --
        # The switched choice is the door that isn't the initial choice and isn't the one the host opened.
        switched_choice = [door for door in doors if door != contestant_choice and door != host_opens][0]
        if switched_choice == prize_door:
            switch_wins += 1

    # Calculate probabilities
    stay_probability = stay_wins / num_trials
    switch_probability = switch_wins / num_trials

    # --- Print Results ---
    print("--- Results ---")
    print("\nStrategy 1: Always STAY with your initial choice.")
    print(f"Total wins: {stay_wins}")
    print("The final probability equation is Wins / Total Trials:")
    print(f"{stay_wins} / {num_trials} = {stay_probability:.4f} (or {stay_probability:.2%})")

    print("\nStrategy 2: Always SWITCH to the other unopened door.")
    print(f"Total wins: {switch_wins}")
    print("The final probability equation is Wins / Total Trials:")
    print(f"{switch_wins} / {num_trials} = {switch_probability:.4f} (or {switch_probability:.2%})")

    print("\n--------------------------------------------------------------")
    print("Conclusion: Yes, switching doors significantly increases your probability of winning.")
    print("You are about twice as likely to win if you switch.")

# Run the simulation
run_monty_hall_simulation()