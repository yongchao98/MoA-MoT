import sys

def solve_caesar_riddle():
    """
    Calculates the length of the longest possible message for Caesar's love letter.
    """

    # Step 1: Explain the plan based on the problem description.
    print("To find the length of the longest message, we must follow these steps:")
    print("1. Define the encryption rule: Letters A-Z are mapped to Roman numerals of their alphabet position (A=1, B=2, etc.). A space remains a space.")
    print("2. Determine the 'cost' of each character, which is the length of its encrypted form.")
    print("3. To maximize the message length, we must use characters with the minimum possible cost.")
    print("4. The maximum message length is the paper's character limit divided by this minimum cost.\n")

    def int_to_roman(num):
        """Converts an integer (1-26) to a Roman numeral string."""
        val = [10, 9, 5, 4, 1]
        syb = ["X", "IX", "V", "IV", "I"]
        roman_num = ''
        i = 0
        # This simplified converter is sufficient for numbers up to 39
        while num > 0:
            for _ in range(num // val[i]):
                roman_num += syb[i]
                num -= val[i]
            i += 1
        return roman_num

    # Step 2: Calculate the cost for each character.
    # The cost for a space is 1.
    min_cost = 1
    min_cost_chars = ['space']
    
    print("Calculating the encryption cost for each letter:")
    for i in range(26):
        letter = chr(ord('A') + i)
        number = i + 1
        roman_representation = int_to_roman(number)
        cost = len(roman_representation)
        # print(f"  - Letter '{letter}' ({number}) -> '{roman_representation}' (Cost: {cost})") # Uncomment for detailed costs
        if cost < min_cost:
            min_cost = cost
            min_cost_chars = [letter]
        elif cost == min_cost:
            # We only need to add the letter if it's not already there (which it won't be)
            min_cost_chars.append(letter)
    
    print("\nCalculation complete.")
    
    # Step 3 & 4: Find the minimum cost and calculate the maximum length.
    paper_capacity = 10000

    print(f"\nThe characters with the lowest encryption cost are: {', '.join(min_cost_chars)}.")
    print(f"The minimum cost to encrypt a single character is: {min_cost}")
    print(f"The total capacity of the paper is: {paper_capacity} characters.")
    
    max_message_length = paper_capacity // min_cost

    print("\nThe final equation to find the longest message length is:")
    print(f"Max Length = (Total Capacity) / (Minimum Cost)")
    print(f"Max Length = {paper_capacity} / {min_cost}")
    print(f"The length of his longest message is {max_message_length}.")

solve_caesar_riddle()
<<<10000>>>