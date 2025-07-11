import sys

def solve():
    """
    Calculates the length of the longest possible message for Caesar.
    """

    # Step 1: Implement a function to convert an integer to a Roman numeral.
    def integer_to_roman(num):
        """Converts an integer to a Roman numeral string."""
        val_map = [
            (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
            (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
            (10, "X"), (9, "IX"), (5, "V"), (4, "IV"),
            (1, "I")
        ]
        roman_numeral = []
        for val, symbol in val_map:
            while num >= val:
                roman_numeral.append(symbol)
                num -= val
        return "".join(roman_numeral)

    # Step 2: Find the minimum possible length for an encrypted character.
    # We assume A=1, B=2, ..., Z=26. We need to find which number in this
    # range produces the shortest Roman numeral.
    min_encrypted_length = sys.maxsize
    
    # We only need to check the numbers corresponding to the alphabet (1-26).
    for i in range(1, 27):
        roman_representation = integer_to_roman(i)
        if len(roman_representation) < min_encrypted_length:
            min_encrypted_length = len(roman_representation)
            
    # Step 3: Define the total capacity of the paper.
    paper_capacity = 10000

    # Step 4: Calculate the maximum message length by dividing the total capacity
    # by the minimum cost per character.
    max_message_length = paper_capacity // min_encrypted_length

    # Step 5: Print the explanation and the final equation.
    print("To write the longest message, Caesar must use characters that have the shortest Roman numeral encryption.")
    print(f"The minimum length of an encrypted character is: {min_encrypted_length}")
    print(f"The total character capacity of the paper is: {paper_capacity}")
    print("\nThe final calculation is:")
    print(f"Paper Capacity / Minimum Encrypted Length = Maximum Message Length")
    print(f"{paper_capacity} / {min_encrypted_length} = {max_message_length}")

solve()
<<<10000>>>