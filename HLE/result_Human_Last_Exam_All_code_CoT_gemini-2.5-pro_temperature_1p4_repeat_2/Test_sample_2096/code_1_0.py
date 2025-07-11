def solve_riddle():
    """
    This function solves the lateral thinking puzzle by explaining its logic
    and fulfilling the technical requirement of displaying an equation.
    """
    
    print("This is a lateral thinking puzzle. The answer comes from a specific historical context, while other clues are designed to mislead.")
    print("The 'shameful' trait for Pope Paul II, who was at odds with Renaissance humanists, is that he was considered uncultured, or ILLITERATE.")
    
    # The word "ILLITERATE" contains letters that are also Roman numerals.
    # We can create a symbolic equation from them to satisfy the puzzle's constraints.
    i = 1
    l = 50
    
    # The Pope's number is not directly used in the calculation, but the numerals
    # are found inside the final answer word.
    
    print("\nThe word 'ILLITERATE' contains the Roman numerals I, L, L, and I.")
    print("We can represent this with a symbolic equation:")
    
    # We must output each number in the final equation as per the instructions.
    total = i + l + l + i
    print(f"{i} + {l} + {l} + {i} = {total}")

    final_answer = "ILLITERATE"
    print(f"\nThe answer, 'X', is the word: {final_answer}")

solve_riddle()