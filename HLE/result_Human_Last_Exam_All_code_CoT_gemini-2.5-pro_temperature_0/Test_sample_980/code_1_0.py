def calculate_next_number():
    """
    This function calculates the next number in the sequence based on a rule
    derived from the prompt.
    """
    word = "Think"
    vowels = "aeiouAEIOU"
    mod_value = 181

    # Find the vowels in the word
    vowels_in_word = [c for c in word if c in vowels]
    
    # Calculate the sum of the ASCII values of the vowels
    vowel_sum = sum(ord(c) for c in vowels_in_word)
    
    # Apply the formula
    result = vowel_sum % mod_value

    # Print the equation step-by-step as requested
    vowel_chars_str = " + ".join([f"ord('{c}')" for c in vowels_in_word])
    vowel_values_str = " + ".join([str(ord(c)) for c in vowels_in_word])
    
    print(f"The next word in the sequence is '{word}'.")
    print("The rule is: (sum of ASCII values of vowels) % 181")
    print(f"Equation: ({vowel_chars_str}) % {mod_value}")
    print(f"= ({vowel_values_str}) % {mod_value}")
    print(f"= {vowel_sum} % {mod_value}")
    print(f"= {result}")

calculate_next_number()