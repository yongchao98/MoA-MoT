def solve_sequence():
    """
    Solves the number sequence puzzle by finding its underlying pattern.
    The sequence is 2, 11, 23, 51, 119, ...
    """
    sequence = [2, 11, 23, 51, 119]
    print(f"Original sequence: {sequence}\n")

    # Step 1: Hypothesize a pattern of the form x_n = 2^(n+1) + k_n, starting from n=2.
    print("Let's test the hypothesis that for n >= 2, the pattern is x_n = 2^(n+1) + k_n.")
    print("We calculate the remainder sequence 'k' for the known terms.")
    
    k_sequence = []
    # We treat the sequence as 1-indexed for the formula, but access it as 0-indexed.
    # The pattern starts from the 2nd term (index 1).
    for i in range(1, len(sequence)):
        n = i + 1 # n-th term in the sequence (n=2, 3, 4, 5)
        power_of_2 = 2**(n+1)
        remainder = sequence[i] - power_of_2
        k_sequence.append(remainder)
        print(f"For x_{n} = {sequence[i]}: k_{n} = {sequence[i]} - 2^({n}+1) = {sequence[i]} - {power_of_2} = {remainder}")

    print(f"\nThe remainder sequence 'k' is: {k_sequence}\n")

    # Step 2: Find the pattern in the remainder sequence 'k'.
    print("Now, let's find the pattern in the remainder sequence by looking at its differences.")
    k_diffs = []
    for i in range(len(k_sequence) - 1):
        diff = k_sequence[i+1] - k_sequence[i]
        k_diffs.append(diff)
        print(f"Difference {i+1}: {k_sequence[i+1]} - {k_sequence[i]} = {diff}")

    print(f"\nThe differences are: {k_diffs}. We can see a geometric progression where each term is multiplied by 3.\n")

    # Step 3: Extrapolate the pattern to find the next term.
    print("The next difference should be the last difference multiplied by 3.")
    last_diff = k_diffs[-1]
    next_diff = last_diff * 3
    print(f"Next difference = {last_diff} * 3 = {next_diff}")

    print("\nNow we can find the next term in the 'k' sequence.")
    last_k = k_sequence[-1]
    next_k = last_k + next_diff
    print(f"Next k = {last_k} + {next_diff} = {next_k}")

    # Step 4: Calculate the final answer.
    # The next term is the 6th term (n=6).
    n_next = len(sequence) + 1
    next_power_of_2 = 2**(n_next + 1)
    next_term = next_power_of_2 + next_k

    print("\nFinally, we can calculate the next term in the original sequence (x_6).")
    print(f"The formula is x_6 = 2^(6+1) + k_6 = 2^7 + {next_k}")
    print(f"So, the final calculation is: {next_power_of_2} + {next_k} = {next_term}")

solve_sequence()
<<<A>>>