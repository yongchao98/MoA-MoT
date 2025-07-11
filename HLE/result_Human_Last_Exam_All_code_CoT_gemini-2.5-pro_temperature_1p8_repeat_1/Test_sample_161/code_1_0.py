def int_to_roman(num):
    """Converts an integer to a Roman numeral string."""
    val_map = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
        (1, "I")
    ]
    roman_str = []
    for val, numeral in val_map:
        # Greedily append the largest possible numeral
        count = num // val
        roman_str.append(numeral * count)
        num -= val * count
    return "".join(roman_str)

def solve_longest_message():
    """
    Calculates the length of the longest possible message.
    """
    # The message is written on paper that can hold 10000 characters.
    paper_capacity = 10000

    # We need to find the minimum possible length of an encrypted character.
    # We check the lengths for the letters A-Z (represented by numbers 1-26).
    min_encrypted_length = float('inf')
    
    # In ASCII, 'A' is 65. We check numbers 1 through 26.
    for i in range(1, 27):
        roman_numeral = int_to_roman(i)
        length = len(roman_numeral)
        if length < min_encrypted_length:
            min_encrypted_length = length
            
    # The minimum length for an encrypted character is the one found above.
    # For example, A=1 -> "I" (length 1), E=5 -> "V" (length 1).
    # To get the longest message, we assume every character in the original message
    # (including spaces) can be encoded using this minimum length.
    
    # Calculate the longest possible message length.
    longest_message = paper_capacity // min_encrypted_length

    # Print the equation as requested.
    print(f"The total capacity of the paper is {paper_capacity} characters.")
    print(f"The minimum cost to encrypt a single character (e.g., 'A' -> 'I') is {min_encrypted_length} character.")
    print("Therefore, the maximum length of the original message is calculated as:")
    print(f"{paper_capacity} / {min_encrypted_length} = {longest_message}")

solve_longest_message()