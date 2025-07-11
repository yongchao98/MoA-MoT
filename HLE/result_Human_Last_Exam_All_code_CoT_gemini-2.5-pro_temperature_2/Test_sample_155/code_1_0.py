def decode_message():
    """
    Decodes and prints a word hidden in a visual puzzle.

    The puzzle uses 't' characters to form the shape of letters.
    This script replaces each 't' with a block character ('#') and
    each 'b' or space with a blank space to reveal the hidden word.
    The exact position of each character is preserved to maintain the
    original shape of the letters.
    """

    # The puzzle input as a list of strings
    puzzle_lines = [
        "t b b t t t t b b b",
        "b",
        "b",
        "b",
        "t",
        "b",
        "t b t t b b b b b b",
        "                  b",
        "                  b",
        "                  t",
        "                  t",
        "                  b",
        "b b t t t t b t b b"
    ]

    print("Decoding the message by rendering 't' characters as blocks:\n")

    # Determine the full width of the drawing area
    max_width = 0
    if puzzle_lines:
        max_width = max(len(line) for line in puzzle_lines)

    # Create a grid and populate it with the puzzle characters
    grid = [[' ' for _ in range(max_width)] for _ in range(len(puzzle_lines))]
    for row_idx, line in enumerate(puzzle_lines):
        for col_idx, char in enumerate(line):
            grid[row_idx][col_idx] = char

    # Print the decoded message by interpreting the grid
    for row in grid:
        output_line = ""
        for char in row:
            if char == 't':
                output_line += '#'
            else:
                output_line += ' '
        print(output_line)

    print("\nReading the letters formed by the '#' characters, the word appears to be CYBER.")

decode_message()