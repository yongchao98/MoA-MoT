import math

def solve_caesar_cipher():
    """
    Calculates the length of the longest possible message under Caesar's
    Roman numeral encryption scheme.
    """

    def int_to_roman(num):
        """Converts an integer to a Roman numeral string."""
        # This is a standard algorithm for integer to Roman numeral conversion.
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]
        roman_num = ""
        for val, numeral in val_map:
            roman_num += numeral * (num // val)
            num %= val
        return roman_num

    # Step 1: Find the shortest Roman numeral representation for any character.
    # We check numbers 1-27 (A-Z and space).
    min_len = float('inf')
    for i in range(1, 28):
        roman_len = len(int_to_roman(i))
        if roman_len < min_len:
            min_len = roman_len
    
    print(f"Step 1: The shortest Roman numeral representation for a single character (e.g., for 'A' mapped to 1 -> 'I') has a length of {min_len}.")

    # Step 2: Define variables and set up the inequality.
    # The message is of length 'n'.
    # A delimiter of length 1 is assumed to separate each Roman numeral.
    paper_capacity = 10000
    delimiter_len = 1
    
    print(f"\nStep 2: Let the message length be 'n'. With a delimiter of length {delimiter_len}, the total encrypted length must be at most {paper_capacity}.")
    print("The governing inequality is:")
    print(f"n * (shortest Roman length) + (n - 1) * (delimiter length) <= {paper_capacity}")
    # This simplifies to n * (min_len + delimiter_len) - delimiter_len <= paper_capacity
    total_cost_per_char = min_len + delimiter_len

    print(f"n * {min_len} + (n - 1) * {delimiter_len} <= {paper_capacity}")
    
    # Step 3: Solve the inequality for 'n'.
    print("\nStep 3: Solving the inequality for 'n':")
    print(f"{min_len}*n + {delimiter_len}*n - {delimiter_len} <= {paper_capacity}")
    print(f"{total_cost_per_char}*n - {delimiter_len} <= {paper_capacity}")
    
    numerator = paper_capacity + delimiter_len
    print(f"{total_cost_per_char}*n <= {paper_capacity} + {delimiter_len}")
    print(f"{total_cost_per_char}*n <= {numerator}")
    
    max_len_float = numerator / total_cost_per_char
    print(f"n <= {numerator} / {total_cost_per_char}")
    print(f"n <= {max_len_float}")
    
    # The length 'n' must be an integer.
    max_len_int = math.floor(max_len_float)
    
    print(f"\nAs the message length must be a whole number, the longest possible message is {max_len_int} characters long.")

solve_caesar_cipher()
<<<5000>>>