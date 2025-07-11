import math
from collections import Counter
import sys

def solve_tverberg_partitions():
    """
    This function calculates the number of ways to partition a set of n items
    into unlabeled subsets of specified sizes. This corresponds to the number
    of Tverberg partitions under the conditions given in the problem.

    The formula is: N = (1 / (k1! * k2! * ...)) * (n! / (a1! * a2! * ...))
    where:
    - n is the total number of items (points).
    - a_i are the sizes of the parts in the partition.
    - k_i are the counts of identical part sizes.
    """

    # For user interaction, the script prompts for input.
    # If run in a non-interactive environment, it uses a default example.
    if sys.stdin.isatty():
        print("Please enter the sizes of the partition, separated by spaces (e.g., '3 3 2 1 1'):")
    input_str = sys.stdin.readline().strip()

    # Use a default example if no input is provided.
    if not input_str:
        print("No input provided. Using example partition '3 3 2 1 1'.")
        input_str = "3 3 2 1 1"

    try:
        parts = [int(p) for p in input_str.split()]
        if not parts:
            print("Error: The partition must have at least one part.", file=sys.stderr)
            return
        if any(p <= 0 for p in parts):
            print("Error: Part sizes must be positive integers.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input. Please enter space-separated integers.", file=sys.stderr)
        return

    # Calculate n = total number of points.
    n = sum(parts)

    # Calculate the multinomial coefficient for labeled groups.
    try:
        n_factorial = math.factorial(n)
        denom_parts_factorial = 1
        for p in parts:
            denom_parts_factorial *= math.factorial(p)
        
        # This will always be an integer division.
        multinomial_coeff = n_factorial // denom_parts_factorial
    except (ValueError, OverflowError):
        print("Error: Calculation resulted in a number too large to handle. Please try a smaller partition.", file=sys.stderr)
        return

    # Count occurrences of each part size to find the k_i values for the correction factor.
    counts = Counter(parts)
    
    # Calculate the correction factor for unlabeled (indistinguishable) groups.
    correction_denom = 1
    for count in counts.values():
        correction_denom *= math.factorial(count)

    # The final result is the multinomial coefficient divided by the correction factor.
    result = multinomial_coeff // correction_denom

    # --- Construct the output string showing the full equation as requested ---

    # Create the multinomial coefficient part of the string, e.g., "10! / (3! * 3! * 2! * 1! * 1!)"
    multinomial_str = f"{n}! / ({' * '.join(f'{p}!' for p in sorted(parts, reverse=True))})"

    # Create the correction factor part of the string, e.g., "1 / (2! * 2!)"
    # We only include counts > 1 for clarity, as 1! = 1.
    counts_for_str = [c for c in counts.values() if c > 1]
    if counts_for_str:
        correction_str = f"1 / ({' * '.join(f'{c}!' for c in sorted(counts_for_str, reverse=True))})"
        # Combine parts into the final equation string.
        full_equation_str = f"({correction_str}) * ({multinomial_str})"
    else:
        # If all part sizes are unique, no correction factor is needed.
        full_equation_str = multinomial_str

    # Print the final, formatted equation with all the numbers and the final result.
    print("\n" + "="*50)
    print("The number of Tverberg's partitions is given by the equation:")
    print(f"N = {full_equation_str}")
    print(f"N = {result}")
    print("="*50)

# Execute the main function.
solve_tverberg_partitions()