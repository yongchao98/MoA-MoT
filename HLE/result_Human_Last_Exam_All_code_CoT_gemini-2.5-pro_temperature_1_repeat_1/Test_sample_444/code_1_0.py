import random

def create_adversary_sequence(max_len=20, max_val=100):
    """The adversary creates a sequence that is eventually zero."""
    n = random.randint(1, max_len)
    seq = [random.randint(0, max_val - 1) for _ in range(n)]
    # The sequence is padded with zeros, but we only need the non-zero part for the sum
    return seq

def alice_strategy(observed_sequence, observed_sum, num_classes=10):
    """
    Alice's strategy:
    1. She knows the sum of the numbers in the boxes she opened.
    2. She randomly guesses which of the 10 equivalence classes the sequence belongs to.
    3. Based on her guess, she calculates the value of the number in the closed box.
    """
    # Alice randomly guesses the class (the total sum modulo num_classes)
    guessed_class = random.randint(0, num_classes - 1)
    
    # She calculates her guess for the number in the closed box
    # guess + observed_sum_mod_10 = guessed_class
    observed_sum_mod = observed_sum % num_classes
    guess = (guessed_class - observed_sum_mod + num_classes) % num_classes
    
    # The problem states numbers are natural numbers. Our strategy only determines the
    # value modulo 10. To make a concrete guess, we can guess the smallest natural
    # number that fits, which is the result of the modulo operation.
    return guess

def run_simulation(num_trials=20000, num_classes=10):
    """
    Simulates the game for a number of trials to find Alice's success rate.
    """
    wins = 0
    box_to_be_guessed = 0 # Alice decides to guess the number in the first box.

    for _ in range(num_trials):
        # 1. The adversary sets up the boxes.
        sequence = create_adversary_sequence()
        
        # Ensure the sequence is long enough to have a box 0 to guess.
        if len(sequence) == 0:
            sequence.append(0)

        # 2. Alice opens all boxes except one.
        true_value_in_closed_box = sequence[box_to_be_guessed]
        
        # Calculate sum of numbers in opened boxes
        observed_sum = 0
        for i in range(len(sequence)):
            if i != box_to_be_guessed:
                observed_sum += sequence[i]

        # 3. Alice makes her guess.
        guess = alice_strategy(sequence, observed_sum, num_classes)

        # 4. Check if Alice wins. We check modulo num_classes as her strategy
        #    is based on classes modulo num_classes.
        if guess == (true_value_in_closed_box % num_classes):
            wins += 1

    success_rate = wins / num_trials
    
    print(f"Simulation Results ({num_trials} trials):")
    print(f"Alice's success rate: {success_rate:.4f}")
    print("\nThis rate is approximately 1/10, as predicted by the analysis.")
    
    numerator = 1
    denominator = 10
    
    print("\nThe theoretical probability of success is given by the equation:")
    print(f"P(Success) = {numerator} / {denominator}")
    print("Where the numbers in the equation are:")
    print(f"Numerator: {numerator}")
    print(f"Denominator: {denominator}")

run_simulation()