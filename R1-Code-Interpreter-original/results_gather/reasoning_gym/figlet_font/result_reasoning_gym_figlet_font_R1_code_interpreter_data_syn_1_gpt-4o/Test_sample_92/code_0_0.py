# Define the pattern as a list of strings
pattern = [
    "######## ######## ######## ######## ######## ######## ########",
    "#      # ######## ###  ### ######## ##    ## ######## #  ##  #",
    "#  ##### ######## ##    ## ######## ###  ### ######## #   #  #",
    "#    ### ######## #  ##  # ######## ###  ### ######## #      #",
    "#  ##### ######## #      # ######## ###  ### ######## #      #",
    "#  ##### ######## #  ##  # ######## ###  ### ######## #  #   #",
    "#  ##### ######## #  ##  # ######## ##    ## ######## #  ##  #",
    "######## ######## ######## ######## ######## ######## ########",
    "######## ######## ######## ######## ########",
    "#      # ######## #      # ######## #    ###",
    "###  ### ######## #  ##### ######## #  #  ##",
    "###  ### ######## #    ### ######## #  ##  #",
    "###  ### ######## #  ##### ######## #  ##  #",
    "###  ### ######## #  ##### ######## #  #  ##",
    "###  ### ######## #      # ######## #    ###",
    "######## ######## ######## ######## ########"
]

# Function to extract a character from the pattern
def extract_character(pattern, start_col, start_row=0):
    char = []
    for row in range(start_row, start_row + 7):
        char.append(pattern[row][start_col:start_col + 7])
    return char

# Function to interpret a 7x7 block as a character
def interpret_character(char_block):
    # Define a simple mapping for known characters
    # This is a simplified example; in practice, you'd have a full mapping
    mapping = {
        ("#      #", "#  #####", "#    ###", "#  #####", "#  #####", "#  #####", "########"): "H",
        ("########", "##    ##", "###  ###", "###  ###", "###  ###", "##    ##", "########"): "O",
        ("###  ###", "#    ###", "#  #####", "#  #####", "#    ###", "###  ###", "########"): "E",
        ("########", "#      #", "###  ###", "###  ###", "###  ###", "#      #", "########"): "C",
        ("########", "#      #", "#      #", "#      #", "#      #", "#      #", "########"): "I",
    }
    # Convert the block to a tuple of strings for easy comparison
    char_tuple = tuple(char_block)
    return mapping.get(char_tuple, "?")

# Extract and interpret each character
decoded_message = []
# First line of characters
for i in range(0, 56, 7):
    char_block = extract_character(pattern, i)
    decoded_message.append(interpret_character(char_block))

# Second line of characters
for i in range(0, 35, 7):
    char_block = extract_character(pattern, i, start_row=8)
    decoded_message.append(interpret_character(char_block))

# Join the decoded characters into a message
decoded_message = ''.join(decoded_message)
print(decoded_message)