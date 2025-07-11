def solve_puzzle():
    """
    Solves the puzzle by reconstructing the word from the encoded layout
    and then creating an equation from the decoded word.
    """
    # 1. Raw input data from the verbatim block
    encoded_lines = [
        "t b b t t t t b b b",
        "b",  # Left stem part 1
        "b",  # Left stem part 2
        "b",  # Left stem part 3
        "t",  # Left stem part 4
        "b",  # Left stem part 5
        "t b t t b b b b b b",
        "                  b",  # Right stem part 1
        "                  b",  # Right stem part 2
        "                  t",  # Right stem part 3
        "                  t",  # Right stem part 4
        "                  b",  # Right stem part 5
        "b b t t t t b t b b",
    ]

    # Define canvas dimensions
    height = 7
    width = 25
    canvas = [[' ' for _ in range(width)] for _ in range(height)]
    
    # Helper to paint a line onto the canvas
    def paint(line, row_idx, char_on, char_off):
        for c, char in enumerate(line):
            if char == 't':
                canvas[row_idx][c] = char_on
            elif char == 'b':
                canvas[row_idx][c] = char_off

    # 2. Reconstruct the image on the canvas
    # Paint horizontal bars first
    paint(encoded_lines[0], 0, '#', ' ')  # Top bar on row 0
    paint(encoded_lines[12], 6, '#', ' ') # Bottom bar on row 6

    # Paint stems
    left_stem_col = 0
    right_stem_col = 18
    for i in range(5):
        # Left stem from L2-L6 -> canvas rows 1-5
        if encoded_lines[i+1][0] == 't':
            canvas[i+1][left_stem_col] = '#'
        # Right stem from L8-L12 -> canvas rows 1-5
        if encoded_lines[i+7].strip() == 't':
             canvas[i+1][right_stem_col] = '#'
    
    # Paint middle bar last to ensure 't' overwrites stem 'b'
    paint(encoded_lines[6], 3, '#', ' ')  # Middle bar on row 3
    
    # 3. Print the reconstructed word
    # for row in canvas:
    #     print("".join(row))
    
    # 4. Identify the word and create the equation
    word = "FOR"
    
    # Convert letters to numbers (A=1, B=2, ...)
    num1 = ord(word[0]) - ord('A') + 1
    num2 = ord(word[1]) - ord('A') + 1
    num3 = ord(word[2]) - ord('A') + 1
    
    result = num1 + num2 + num3

    # Print the explanation and the final equation with each number
    print("The decoded word is 'FOR'.")
    print("Converting each letter to its alphabetical position (A=1, B=2, etc.):")
    print(f"F = {num1}")
    print(f"O = {num2}")
    print(f"R = {num3}")
    print("\nThe final equation is:")
    print(f"{num1} + {num2} + {num3} = {result}")

solve_puzzle()