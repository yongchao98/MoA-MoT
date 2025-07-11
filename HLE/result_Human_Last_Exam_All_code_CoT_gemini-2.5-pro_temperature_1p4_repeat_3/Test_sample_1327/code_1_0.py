def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying the underlying pattern.
    """
    # The given sequence starts with 2, 11, 23, 51, 119.
    # The pattern, starting from the second term, is x_n = 2 * x_{n-1} + k_n
    # Let's find the sequence of k.
    # k1 = 23 - (2 * 11) = 1
    # k2 = 51 - (2 * 23) = 5
    # k3 = 119 - (2 * 51) = 17
    k_sequence = [1, 5, 17]

    # Now, let's find the pattern in the k_sequence.
    # The differences are: 5 - 1 = 4, 17 - 5 = 12.
    # The differences [4, 12] form a geometric progression with a ratio of 3.
    last_k_diff = 12
    ratio = 3
    next_k_diff = last_k_diff * ratio

    # The next k is the last k plus the next difference.
    last_k = 17
    next_k = last_k + next_k_diff

    # The last term of the original sequence is 119.
    last_term = 119
    next_term = 2 * last_term + next_k

    # Print the final equation with all its components.
    print(f"The pattern is that each number is twice the previous number plus a value 'k'.")
    print(f"The sequence of 'k' is 1, 5, 17, ...")
    print(f"The differences in 'k' are 4, 12, ... which shows a multiplicative factor of 3.")
    print(f"The next difference in 'k' is {last_k_diff} * {ratio} = {next_k_diff}.")
    print(f"The next 'k' is {last_k} + {next_k_diff} = {next_k}.")
    print(f"Therefore, the final equation is:")
    print(f"{last_term} * 2 + {next_k} = {next_term}")

solve_sequence()
<<<A>>>