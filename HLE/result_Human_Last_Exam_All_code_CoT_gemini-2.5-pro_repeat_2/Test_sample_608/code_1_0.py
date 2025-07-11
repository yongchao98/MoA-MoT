import random

def monty_hall_simulation(num_trials):
    """
    Simulates the Monty Hall problem for a given number of trials.
    Returns the number of wins for both staying and switching strategies.
    """
    stay_wins = 0
    switch_wins = 0

    for _ in range(num_trials):
        # Doors are represented by numbers 0, 1, 2
        # Let's place the prize behind a random door
        prize_door = random.randint(0, 2)
        
        # Player makes an initial choice
        player_choice = random.randint(0, 2)
        
        # Host opens a door that is not the player's choice and not the prize door
        opened_by_host = next(door for door in range(3) if door != player_choice and door != prize_door)
        
        # Player decides to switch
        # The switched door is the one that is not the player's original choice and not the one opened by the host
        switched_choice = next(door for door in range(3) if door != player_choice and door != opened_by_host)

        # --- Calculate wins for each strategy ---
        
        # 'Stay' strategy wins if the initial choice was the prize
        if player_choice == prize_door:
            stay_wins += 1
            
        # 'Switch' strategy wins if the switched choice is the prize
        if switched_choice == prize_door:
            switch_wins += 1
            
    return stay_wins, switch_wins

# --- Main execution ---
if __name__ == "__main__":
    num_trials = 10000
    
    stay_wins, switch_wins = monty_hall_simulation(num_trials)
    
    stay_probability = stay_wins / num_trials
    switch_probability = switch_wins / num_trials

    print(f"Simulating the Monty Hall problem {num_trials} times...\n")

    print("--- Strategy 1: Always Stay ---")
    print(f"Wins: {stay_wins} out of {num_trials} trials")
    print(f"Win probability: {stay_wins} / {num_trials} = {stay_probability:.4f} or {stay_probability:.2%}\n")
    
    print("--- Strategy 2: Always Switch ---")
    print(f"Wins: {switch_wins} out of {num_trials} trials")
    print(f"Win probability: {switch_wins} / {num_trials} = {switch_probability:.4f} or {switch_probability:.2%}\n")

    print("Conclusion: The simulation shows that switching doors approximately doubles your chances of winning.")

    # The final answer is the ratio of switch-win-probability to stay-win-probability
    # In theory, this is (2/3) / (1/3) = 2
    # We will use the simulated value.
    final_answer = switch_probability / stay_probability if stay_probability > 0 else 0
    # Let's just output the conclusion in text form
    final_answer_text = "Yes, changing it increases my probability of winning from approximately 1/3 to 2/3."
    print(f"\nFinal Answer: {final_answer_text}")

    # The actual question is "Does changing it increase my probability of winning?"
    # The answer is simply "Yes".
    # But I should follow the instructions and give a "final answer" block.
    # The most direct answer to "Does it increase probability?" is "Yes".
    # I will output "Yes" as the final answer.
    # <<<Yes>>>