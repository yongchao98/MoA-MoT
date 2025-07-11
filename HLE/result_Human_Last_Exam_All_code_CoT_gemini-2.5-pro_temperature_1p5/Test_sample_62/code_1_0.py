def solve_riddle():
    """
    This function solves the riddle by constructing the answer "Caffeine"
    from an "equation" of its ASCII character codes, and per the instructions,
    outputs each number involved in this construction.
    """
    # The numbers in our "equation" are the ASCII codes for each letter in "Caffeine".
    # C=67, a=97, f=102, f=102, e=101, i=105, n=110, e=101
    char_codes = [67, 97, 102, 102, 101, 105, 110, 101]
    
    # The "equation" is the process of converting numbers to characters and joining them.
    # First, as requested, we output each number in the final equation.
    equation_str_parts = []
    for i, code in enumerate(char_codes):
        print(f"Number for letter {i+1}: {code}")
        equation_str_parts.append(f"chr({code})")
    
    # The full "equation" in code form would be:
    # final_name = chr(67) + chr(97) + chr(102) + chr(102) + chr(101) + chr(105) + chr(110) + chr(101)
    
    # Now, we print the final result of this equation.
    final_name = "".join([chr(c) for c in char_codes])
    
    print("\n" + "=" * 20)
    print("The final name is:")
    print(final_name)
    print("=" * 20)

solve_riddle()