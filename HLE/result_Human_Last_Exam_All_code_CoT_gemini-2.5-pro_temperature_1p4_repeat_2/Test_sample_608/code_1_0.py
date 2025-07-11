import random

def run_monty_hall_simulation(num_simulations):
    """
    Simulates the Monty Hall problem for a given number of iterations.
    """
    stay_wins = 0
    switch_wins = 0
    
    # The doors are represented by numbers 1, 2, and 3
    doors = [1, 2, 3]

    for _ in range(num_simulations):
        # 1. The prize is placed behind a random door
        prize_door = random.choice(doors)
        
        # 2. The contestant makes an initial choice
        player_choice = random.choice(doors)

        # 3. The host opens a door that is NOT the player's choice and NOT the prize door
        # Create a list of doors the host can open
        host_options = [d for d in doors if d != player_choice and d != prize_door]
        host_opened_door = random.choice(host_options)

        # ---- Evaluate the "Stay" Strategy ----
        # The player wins if their initial choice was the prize door
        if player_choice == prize_door:
            stay_wins += 1

        # ---- Evaluate the "Switch" Strategy ----
        # The player switches to the only other door that is not their original choice
        # and not the one the host opened.
        switched_choice = [d for d in doors if d != player_choice and d != host_opened_door][0]
        if switched_choice == prize_door:
            switch_wins += 1
            
    # --- Print Results ---
    print(f"Running {num_simulations} simulations...\n")

    # Results for staying
    print("--- Strategy: Always Stay ---")
    stay_probability = (stay_wins / num_simulations) * 100
    print(f"Wins: {stay_wins}")
    print(f"Total Games: {num_simulations}")
    print(f"Final Equation: {stay_wins} / {num_simulations} = {stay_probability:.2f}% Probability of Winning")
    print("-" * 35)
    
    # Results for switching
    print("--- Strategy: Always Switch ---")
    switch_probability = (switch_wins / num_simulations) * 100
    print(f"Wins: {switch_wins}")
    print(f"Total Games: {num_simulations}")
    print(f"Final Equation: {switch_wins} / {num_simulations} = {switch_probability:.2f}% Probability of Winning")
    print("-" * 35)

    print("\nConclusion: The simulation shows that switching doors results in winning approximately 2/3 of the time, while staying only wins about 1/3 of the time.")


# Run the simulation with 10,000 games
if __name__ == "__main__":
    run_monty_hall_simulation(10000)
<<<Yes>>>