import math

def main():
    """
    Calculates the length of the longest possible message Caesar can write.
    """

    # 1. Define character-to-number mapping and find the shortest Roman numeral representation.
    # We assume A-Z maps to 1-26, and space maps to 27.
    # To find the character with the shortest encrypted form, we need to find the
    # length of the Roman numeral for each number from 1 to 27.

    # A helper function to calculate the length of a Roman numeral for a given integer.
    def get_roman_numeral_length(n):
        """Calculates the length of the Roman numeral representation of n."""
        if not 0 < n < 4000:
            # The standard Roman numeral system does not represent 0 or numbers >= 4000.
            # Numbers 1-27 are well within this range.
            raise ValueError("Input must be between 1 and 3999")
        
        # Pairs of (value, length of Roman symbol)
        # E.g., 900 is "CM" (length 2), 500 is "D" (length 1)
        val_map = [
            (1000, 1), (900, 2), (500, 1), (400, 2), (100, 1), (90, 2),
            (50, 1), (40, 2), (10, 1), (9, 2), (5, 1), (4, 2), (1, 1)
        ]
        
        length = 0
        for val, l in val_map:
            count = n // val
            length += count * l
            n %= val
        return length

    # Find the minimum length of an encrypted character.
    min_len = float('inf')
    for i in range(1, 28): # For numbers 1-27
        length = get_roman_numeral_length(i)
        if length < min_len:
            min_len = length

    print(f"To write the longest message, Caesar should use characters that have the shortest encrypted form.")
    print(f"The shortest Roman numeral representation for any character (A-Z, space) has a length of: {min_len}\n")

    # 2. Formulate and solve the equation for the maximum message length (N).
    paper_capacity = 10000
    separator_len = 1

    print("The total length of the encrypted message is determined by the sum of the lengths of N Roman numerals")
    print("and the N-1 separators between them.\n")
    print("Formula: Total Length = N * (length of each numeral) + (N - 1) * (length of each separator)")
    print("Using the shortest possible numeral and a single-character separator:\n")

    print("Total Length = N * min_len + (N - 1) * separator_len <= paper_capacity")
    print("Substituting the values:")
    print(f"N * {min_len} + (N - 1) * {separator_len} <= {paper_capacity}")
    
    # Simplify the inequality
    print("\nSimplifying the inequality:")
    min_len_plus_sep = min_len + separator_len
    paper_capacity_plus_sep = paper_capacity + separator_len
    
    print(f"N * {min_len} + N - {separator_len} <= {paper_capacity}")
    print(f"N * ({min_len} + {separator_len}) <= {paper_capacity} + {separator_len}")
    print(f"N * {min_len_plus_sep} <= {paper_capacity_plus_sep}")

    # Solve for N
    max_n_float = paper_capacity_plus_sep / min_len_plus_sep
    max_n_int = math.floor(max_n_float)

    print("\nSolving for N:")
    print(f"N <= {paper_capacity_plus_sep} / {min_len_plus_sep}")
    print(f"N <= {max_n_float}")
    
    print("\nSince the length of the message must be an integer, we take the floor of the result.")
    print(f"The maximum length of the message is {max_n_int}.")

if __name__ == "__main__":
    main()
<<<5000>>>