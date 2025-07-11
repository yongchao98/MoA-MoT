import random

def alice_strategy_simulation():
    """
    Simulates Alice's strategy and its failure against a worst-case sequence.
    """
    # 1. SETUP
    # Let's define a hypothetical representative sequence `r`.
    # For this example, let's say it's the sequence of even numbers.
    # In reality, this is one of uncountably many representatives from the set R.
    r = [2 * i for i in range(100)] # We only need a prefix

    # Alice decides her strategy: leave the first M boxes closed.
    M = 10
    closed_box_indices = list(range(M))
    print(f"Alice's Strategy: Leave boxes 0 through {M-1} closed.")
    print("-" * 20)

    # 2. ADVERSARY'S MOVE
    # The adversary, knowing Alice's strategy, creates a malicious sequence `s`.
    # It's based on `r` but differs on all of Alice's chosen closed boxes.
    s = list(r) # Make a copy
    for i in closed_box_indices:
        s[i] = r[i] + 1 # Make it different
    
    print("A malicious sequence `s` is created.")
    print(f"Representative `r` starts with: {r[:M]}")
    print(f"Malicious sequence `s` starts with: {s[:M]}")
    print("-" * 20)

    # 3. ALICE'S TURN
    # Alice opens the boxes *not* in her closed set.
    # From these, she correctly identifies the representative `r`.
    # (We simulate this by just giving her `r`).
    identified_r = r
    print("Alice opens the other boxes and correctly identifies the representative `r`.")

    # Alice makes her guess. She picks a random closed box.
    box_to_guess_index = random.choice(closed_box_indices)
    print(f"Alice randomly chose to guess the number in box {box_to_guess_index}.")

    # Her guess is the value from the representative sequence.
    alices_guess = identified_r[box_to_guess_index]
    
    # We check the actual value in the box.
    actual_value = s[box_to_guess_index]

    print(f"Alice's guess for box {box_to_guess_index}: {alices_guess}")
    print(f"Actual value in box {box_to_guess_index}: {actual_value}")
    print("-" * 20)

    # 4. OUTCOME
    if alices_guess == actual_value:
        print("Result: Alice's guess was CORRECT.")
        success = True
    else:
        print("Result: Alice's guess was WRONG.")
        success = False
        
    # For this malicious sequence, her probability of success is 0.
    # We can show the simple equation for her success chance with this sequence.
    num_closed_boxes = len(closed_box_indices)
    num_correct_guesses_possible = 0
    for i in closed_box_indices:
        if identified_r[i] == s[i]:
            num_correct_guesses_possible += 1
            
    print("\nAnalysis of Alice's chances for this specific sequence:")
    print(f"Number of winning choices: {num_correct_guesses_possible}")
    print(f"Total number of choices: {num_closed_boxes}")
    print(f"Success Probability = {num_correct_guesses_possible} / {num_closed_boxes} = {num_correct_guesses_possible/num_closed_boxes}")

if __name__ == '__main__':
    alice_strategy_simulation()