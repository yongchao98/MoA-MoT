def solve_sequence():
    """
    This function solves the number sequence puzzle.
    """
    # The rule is that each number is calculated from a word in the prompt's text.
    # After analysis, the proposed formula is:
    # Value = (Number of Vowels) * 31 + (Number of Consonants)
    
    # The next word in the sequence is "the".
    word = "the"
    
    vowels = "aeiouAEIOU"
    num_vowels = 0
    num_consonants = 0
    
    for char in word:
        if char in vowels:
            num_vowels += 1
        elif char.isalpha():
            num_consonants += 1
            
    # These are K and R from the formula: Value = K * 31 + R
    K = num_vowels
    R = num_consonants
    
    # Calculate the final number
    next_number = K * 31 + R
    
    print(f"The next word in the source text is '{word}'.")
    print(f"The number of vowels (K) is: {K}")
    print(f"The number of consonants (R) is: {R}")
    print("The formula is: Value = K * 31 + R")
    print(f"The final equation is: {next_number} = {K} * 31 + {R}")

solve_sequence()