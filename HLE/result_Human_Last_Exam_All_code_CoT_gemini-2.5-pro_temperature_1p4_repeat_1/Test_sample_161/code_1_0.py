def main():
    """
    This script calculates the length of the longest possible message Caesar can write
    under the given constraints.
    """

    # The total number of characters the paper can hold.
    total_capacity = 10000

    def int_to_roman(num: int) -> str:
        """Converts an integer to its standard Roman numeral representation."""
        val_map = [
            (50, 'L'), (40, 'XL'), (10, 'X'), (9, 'IX'),
            (5, 'V'), (4, 'IV'), (1, 'I')
        ]
        roman_numeral = ""
        # The relevant characters A-Z and space map to numbers 1-27.
        # We only need to handle numbers in this range.
        for val, numeral in val_map:
            while num >= val:
                roman_numeral += numeral
                num -= val
        return roman_numeral

    # We assume A-Z map to 1-26 and space maps to 27.
    # To find the longest message, we need the shortest possible encrypted character.
    # We will find the minimum length by checking all possible characters.
    min_encrypted_length = float('inf')

    # Iterate through numbers 1 to 27 to find the shortest Roman numeral.
    for i in range(1, 28):
        roman = int_to_roman(i)
        if len(roman) < min_encrypted_length:
            min_encrypted_length = len(roman)

    # Calculate the maximum possible message length.
    # We use integer division as we can't write a fraction of a character.
    max_message_length = total_capacity // min_encrypted_length
    
    print("To write the longest message, we must use characters that have the shortest encrypted length.")
    print(f"The minimum length of an encrypted character (e.g., 'A' -> 'I', 'E' -> 'V', 'J' -> 'X') is {min_encrypted_length}.")
    print("\nThe paper's total character capacity is 10000.")
    print("The length of the longest message is calculated by dividing the total capacity by the minimum cost per character.")
    print("\nFinal Equation:")
    print(f"{total_capacity} / {min_encrypted_length} = {max_message_length}")

if __name__ == "__main__":
    main()