def solve_caesar_cipher_length():
    """
    Calculates the length of the longest possible message Caesar can write.
    """
    # The paper can hold a maximum of 10000 characters for the encrypted message.
    paper_capacity = 10000

    # We assume each letter A-Z is converted to its position 1-26,
    # and then to a Roman numeral.
    # A dictionary mapping numbers 1-26 to their Roman numeral representation.
    roman_map = {
        1: "I", 2: "II", 3: "III", 4: "IV", 5: "V",
        6: "VI", 7: "VII", 8: "VIII", 9: "IX", 10: "X",
        11: "XI", 12: "XII", 13: "XIII", 14: "XIV", 15: "XV",
        16: "XVI", 17: "XVII", 18: "XVIII", 19: "XIX", 20: "XX",
        21: "XXI", 22: "XXII", 23: "XXIII", 24: "XXIV", 25: "XXV",
        26: "XXVI"
    }

    # Calculate the encryption cost (length) for each letter.
    letter_costs = [len(roman_map[i]) for i in range(1, 27)]

    # The cost for a space character is assumed to be 1 (the most efficient).
    space_cost = 1

    # To find the longest message, we must use characters with the minimum encryption cost.
    min_letter_cost = min(letter_costs)
    min_cost_per_char = min(min_letter_cost, space_cost)

    # The longest message is composed entirely of characters with this minimum cost.
    # Its length is the total capacity divided by the minimum cost.
    max_message_length = paper_capacity // min_cost_per_char

    print("To find the longest message, we must use characters that have the shortest encrypted form.")
    print(f"The total capacity of the paper is {paper_capacity} characters.")
    print(f"The minimum cost to encrypt a single character (e.g., 'A' -> 'I') is {min_cost_per_char} character.")
    print("\nThe length of the longest message is calculated as:")
    print(f"{paper_capacity} (Total Capacity) / {min_cost_per_char} (Minimum Cost) = {max_message_length}")
    print(f"\nTherefore, the length of his longest message is {max_message_length} characters.")

solve_caesar_cipher_length()
<<<10000>>>