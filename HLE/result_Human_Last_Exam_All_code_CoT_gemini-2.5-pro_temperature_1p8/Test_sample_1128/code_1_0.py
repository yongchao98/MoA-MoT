import math

def count_cool_strings():
    """
    Calculates the number of "cool strings" of maximal length 3n.

    The logic follows a proof by contradiction:
    1. A cool string of length 3n must have each of the n symbols appear exactly 3 times.
    2. The structure of such a string must be P1*P2*P3, where P1, P2, P3 are permutations
       containing the first, second, and third occurrences of the symbols, respectively.
    3. This structure is only cool if for every symbol S, S is the last symbol of P1
       and the first symbol of P3.
    4. This is only possible if n=1. For n>1, it's a contradiction, so no such strings exist.
    """
    try:
        n_input = input("Enter the number of symbols (n): ")
        n = int(n_input)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a positive integer.")
        return

    print(f"For n = {n}:")
    if n == 1:
        # For n=1, there is one symbol 'A'. The only string of length 3 is 'AAA'.
        # Between the 1st and 3rd 'A' is one 'A'. The sequence of symbols is ('A'). No duplicates.
        # This string is cool. Number of strings is 1.
        result = 1
        equation_str = "1"
        print("The number of cool strings of maximal length is 1.")
    else:
        # For n > 1, the logical contradiction holds. No such string exists.
        result = 0
        equation_str = "0"
        print("The number of cool strings of maximal length is 0.")

    print("\nThe final equation is:")
    # Per instructions, print each number of the final equation.
    # The final equation is simply the result itself.
    for char in equation_str:
        if char.isdigit():
            print(char)

# Execute the function
count_cool_strings()