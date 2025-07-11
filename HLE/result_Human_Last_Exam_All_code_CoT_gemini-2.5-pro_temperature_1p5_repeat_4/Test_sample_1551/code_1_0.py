def solve_puzzle():
    """
    This function calculates and prints the answers to the user's questions based on mathematical derivations.
    """

    # --- Question 1: What is the expected value of X_19? ---
    n_19 = 19
    # For odd n, the formula for the expected value is (n-1)/2.
    expected_value_19 = (n_19 - 1) / 2
    print(f"1. The expected value of X_19 is calculated using the formula for odd n: (19 - 1) / 2 = {int(expected_value_19)}.")

    # --- Question 2: What is the expected value of X_20? ---
    # For even n, the parity of the distance prevents the end condition from ever being met.
    print("2. The expected value of X_20 is infinite, as the game never ends for an even number of friends. We can represent this as inf.")

    # --- Question 3: In general, for odd n>1 what is the expected value of X_n? ---
    print("3. For a general odd n > 1, the expected value of X_n is given by the formula: (n - 1) / 2.")

    # --- Question 4: For odd n > 30 what is the expected number of times that between the two gifts there are exactly 10 and n-12 friends? ---
    # This corresponds to the unoriented distance being 11. The expected number of visits is 2.
    expected_visits = 2
    print(f"4. For any odd n > 30, the expected number of times the specified state occurs is {expected_visits}.")
    
    # --- Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value ? ---
    # Yes, for odd n, the random walk of the distance is an irreducible Markov chain on a finite state space,
    # with absorption states that are reachable from any other state.
    answer_q5 = "Yes"
    print(f"5. It is true that the game ends in a finite amount of time with probability 1 for every odd n > 1. The answer is: {answer_q5}.")

solve_puzzle()