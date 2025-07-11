import io

def decode_word():
    """
    Decodes a word from a visual pattern of 'b' and 't' characters.
    'b' is treated as a filled block ('#').
    't' is treated as an empty block (' ').
    """
    encoded_text = r'''t b b t t t t b b b
b
b
b
t
b
t b t t b b b b b b
                  b
                  b
                  t
                  t
                  b
b b t t t t b t b b'''

    # Use StringIO to handle multiline string and line endings consistently
    file_like_object = io.StringIO(encoded_text)
    lines = file_like_object.readlines()
    
    # Strip any trailing newline characters from each line
    lines = [line.rstrip('\n') for line in lines]

    # Find the maximum width of the pattern to align all lines
    max_width = 0
    for line in lines:
        if len(line) > max_width:
            max_width = len(line)

    print("--- Decoded Pattern ---")
    
    # Process and print each line to reveal the visual word
    for line in lines:
        # Pad the line with spaces to the maximum width for alignment
        padded_line = line.ljust(max_width)
        
        visual_line = ""
        for char in padded_line:
            if char == 'b':
                visual_line += '#'
            else: # Treat 't' and spaces as empty
                visual_line += ' '
        print(visual_line)
        
    print("-" * max_width)

# Run the decoding function
decode_word()
