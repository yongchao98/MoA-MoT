def solve_caesar_cipher():
    """
    Calculates the length of the longest message Caesar can write using his
    new Roman numeral encryption.
    """

    def int_to_roman(n):
        """Converts an integer (from 1 to 27) to its Roman numeral representation."""
        # This function only needs to handle numbers up to 27, but is built to be general.
        val_map = [
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")
        ]
        
        tens = n // 10
        remainder = n % 10
        
        roman_num = "X" * tens
        
        for val, numeral in val_map:
            if remainder >= val:
                roman_num += numeral
                remainder -= val
                break # This break is for the subtractive cases like 9 and 4
        
        roman_num += "I" * remainder
        
        return roman_num

    # We assume the 27 characters (A-Z and space) are mapped to numbers 1-27.
    # To find the longest message, we must find the shortest possible encrypted character length.
    min_encrypted_length = float('inf')
    
    # Iterate through all possible character values (1 to 27)
    for i in range(1, 28):
        roman_representation = int_to_roman(i)
        current_length = len(roman_representation)
        if current_length < min_encrypted_length:
            min_encrypted_length = current_length

    # The total capacity of the paper is given.
    paper_capacity = 10000

    # The maximum length of the original message is the total capacity divided by the 
    # length of the shortest possible encrypted character.
    max_message_length = paper_capacity // min_encrypted_length

    print("To create the longest message, Caesar must use characters that have the shortest encrypted length.")
    print(f"The minimum possible length for a single encrypted character is: {min_encrypted_length}")
    print(f"The total capacity of the paper is: {paper_capacity} characters")
    print("\nThe final equation to find the maximum message length is:")
    print(f"{paper_capacity} / {min_encrypted_length} = {max_message_length}")

# Execute the solution
solve_caesar_cipher()

<<<10000>>>