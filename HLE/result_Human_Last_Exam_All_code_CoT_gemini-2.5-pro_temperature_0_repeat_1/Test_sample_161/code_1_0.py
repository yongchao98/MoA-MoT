import sys

def number_to_roman(num):
    """Converts a positive integer to its Roman numeral representation."""
    # This function handles numbers from 1 to 3999.
    # For this problem, we only need to handle numbers 1-26.
    val_map = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    
    roman_numeral = ""
    for val, numeral in val_map:
        # Add the Roman numeral symbol for the value as many times as it fits.
        while num >= val:
            roman_numeral += numeral
            num -= val
    return roman_numeral

def solve_caesar_cipher_limit():
    """
    Calculates the maximum length of a message given the encryption rules.
    """
    # Step 1: Find the minimum possible length of a single encrypted character.
    # We assume A=1, B=2, ..., Z=26. The shortest Roman numeral for any of these
    # will be our minimum encoding length.
    min_encoding_length = sys.maxsize 

    # Iterate through numbers 1 to 26, representing letters A to Z.
    for i in range(1, 27):
        roman = number_to_roman(i)
        length = len(roman)
        
        if length < min_encoding_length:
            min_encoding_length = length
    
    # The shortest possible Roman numeral for any positive integer is 1 (e.g., I, V, X).
    # Since we found that A=1 becomes 'I' (length 1), our minimum length is 1.
    # We don't need to know the space's encoding; if it's longer than 1, we just
    # wouldn't use it in the message.

    # Step 2: Define the total character capacity of the paper.
    paper_capacity = 10000

    # Step 3: Calculate the maximum message length.
    max_message_length = paper_capacity // min_encoding_length

    # Step 4: Print the explanation and the final equation with its numbers.
    print("To find the longest possible message, we divide the total paper capacity by the minimum length of a single encrypted character.")
    print(f"Total paper capacity: {paper_capacity}")
    print(f"Minimum character encoding length: {min_encoding_length}")
    print("\nFinal equation:")
    print(f"{paper_capacity} / {min_encoding_length} = {max_message_length}")

solve_caesar_cipher_limit()