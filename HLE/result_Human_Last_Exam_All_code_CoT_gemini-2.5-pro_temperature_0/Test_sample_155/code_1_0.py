def decode_message():
    """
    This function reconstructs a 2D canvas from fragmented text input
    to reveal a hidden word.
    """
    # Define the canvas size. Width is determined by the longest line and indentation.
    width = 30
    height = 13
    canvas = [[' ' for _ in range(width)] for _ in range(height)]

    # The input is provided as fragments with specific locations.
    # Horizontal fragments (row, content)
    horizontal_fragments = [
        (0, "t b b t t t t b b b"),
        (6, "t b t t b b b b b b"),
        (12, "b b t t t t b t b b")
    ]

    # Vertical fragments (column, start_row, content)
    vertical_fragments = [
        (0, 1, "bbtb"), # Note: The 't' is the 4th character in this sequence
        (18, 7, "bbtb") # Note: The two 't's are the 3rd and 4th characters
    ]

    # Draw horizontal fragments onto the canvas
    for row, content in horizontal_fragments:
        for i, char in enumerate(content):
            if i < width:
                canvas[row][i] = char

    # Draw vertical fragments onto the canvas
    # Left vertical segment
    canvas[1][0] = 'b'
    canvas[2][0] = 'b'
    canvas[3][0] = 'b'
    canvas[4][0] = 't'
    canvas[5][0] = 'b'

    # Right vertical segment (with indentation)
    col_index = 18
    canvas[7][col_index] = 'b'
    canvas[8][col_index] = 'b'
    canvas[9][col_index] = 't'
    canvas[10][col_index] = 't'
    canvas[11][col_index] = 'b'

    # Print the reconstructed canvas
    print("Reconstructed message:")
    for row in canvas:
        print("".join(row))

decode_message()