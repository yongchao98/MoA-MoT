def solve_character_riddle():
    """
    This script deciphers a riddle about a Chinese character's strokes.
    It identifies the character and prints an "equation" based on its stroke count.
    """
    # The riddle describes the components of the character.
    # "One horizontal stroke, another horizontal stroke..." implies two horizontal strokes.
    num_horizontal_strokes = 2

    # "...one vertical on the left, one vertical on the right..." implies two vertical strokes.
    num_vertical_strokes = 2

    # The character is formed by combining these strokes.
    # 2 horizontal strokes + 2 vertical strokes = a character with 4 strokes.
    # This character is 井 (jǐng), which means "well".
    total_strokes = num_horizontal_strokes + num_vertical_strokes
    final_character = "井"

    print("Based on the riddle, the character is composed of:")
    print(f"- {num_horizontal_strokes} horizontal strokes")
    print(f"- {num_vertical_strokes} vertical strokes")
    print("\nThe stroke equation for this character is:")
    
    # As requested, printing each number in the final equation.
    print(f"{num_horizontal_strokes} + {num_vertical_strokes} = {total_strokes}")
    
    print(f"\nThe character referred to by the riddle is: {final_character}")

solve_character_riddle()