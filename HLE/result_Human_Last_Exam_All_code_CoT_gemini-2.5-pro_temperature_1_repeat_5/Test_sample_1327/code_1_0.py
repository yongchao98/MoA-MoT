def solve_sequence():
    """
    This function solves the number sequence puzzle by finding the underlying patterns.
    """
    # The given sequence
    s = [2, 11, 23, 51, 119]

    print("Analyzing the sequence: 2, 11, 23, 51, 119, ...")
    print("Let's test the pattern: next_number = 2 * current_number + k")
    print("-" * 30)

    # Calculate the sequence of k
    k = []
    for i in range(len(s) - 1):
        k_val = s[i+1] - 2 * s[i]
        k.append(k_val)
        print(f"{s[i+1]} = 2 * {s[i]} + {k_val}")

    print("-" * 30)
    print(f"The sequence of k is: {k}")
    print("Let's find the pattern in the k sequence [7, 1, 5, 17], ignoring the first term.")
    print("The pattern seems to be: next_k = 3 * current_k + 2")

    # Verify the pattern for k
    print("Verifying the pattern for k:")
    print(f"{k[2]} = 3 * {k[1]} + 2  =>  5 = 3 * 1 + 2. This is correct.")
    print(f"{k[3]} = 3 * {k[2]} + 2  =>  17 = 3 * 5 + 2. This is correct.")
    print("-" * 30)

    # Calculate the next k
    next_k = 3 * k[-1] + 2
    print(f"Now, let's calculate the next k:")
    print(f"next_k = 3 * {k[-1]} + 2 = {next_k}")
    print("-" * 30)

    # Calculate the next number in the main sequence
    last_s = s[-1]
    next_s = 2 * last_s + next_k
    print("Finally, let's find the next number in the main sequence:")
    print(f"next_number = 2 * {last_s} + {next_k} = {next_s}")
    print("-" * 30)
    print(f"The number that should be filled in the parentheses is {next_s}.")

solve_sequence()