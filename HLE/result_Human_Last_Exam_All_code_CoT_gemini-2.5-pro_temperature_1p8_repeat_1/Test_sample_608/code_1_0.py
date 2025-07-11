import random

def run_monty_hall_simulation():
    """
    Simulates the Monty Hall problem to demonstrate the probabilities
    of winning by either staying or switching doors.
    """
    num_simulations = 10000
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_simulations):
        # Doors are represented by numbers 0, 1, 2
        doors = [0, 1, 2]
        
        # Randomly place the prize behind one door
        prize_door = random.choice(doors)
        
        # Player makes an initial random choice
        player_choice = random.choice(doors)

        # The host must open a door that is not the player's choice
        # and not the prize door.
        doors_host_can_open = []
        for door in doors:
            if door != player_choice and door != prize_door:
                doors_host_can_open.append(door)
        
        host_opened_door = random.choice(doors_host_can_open)

        # If the player stays, they win if their initial choice was the prize.
        if player_choice == prize_door:
            stay_wins += 1

        # The player can switch to the remaining closed door. This is the door
        # that is not their initial choice and not the one the host opened.
        # A simpler way to think about the switch is: if the player's initial
        # choice was NOT the prize, switching will ALWAYS win.
        if player_choice != prize_door:
            switch_wins += 1

    # --- Output the results ---
    print(f"Running {num_simulations} simulations of the Monty Hall game...")
    print("-" * 50)
    
    print("Strategy 1: Always STAY with the initial door.")
    print(f"Wins = {stay_wins}")
    print(f"Total Games = {num_simulations}")
    # Final equation for staying
    print(f"Win Probability by Staying: {stay_wins} / {num_simulations} = {stay_wins / num_simulations:.4f}")
    print("Theoretical Probability: 1 / 3 = 0.3333\n")
    
    print("Strategy 2: Always SWITCH to the other door.")
    print(f"Wins = {switch_wins}")
    print(f"Total Games = {num_simulations}")
    # Final equation for switching
    print(f"Win Probability by Switching: {switch_wins} / {num_simulations} = {switch_wins / num_simulations:.4f}")
    print("Theoretical Probability: 2 / 3 = 0.6667\n")

    print("Conclusion: Yes, switching doors doubles your probability of winning.")


if __name__ == '__main__':
    run_monty_hall_simulation()