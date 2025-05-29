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

# Extract and print each character block for manual inspection
for i in range(0, 56, 7):
    char_block = extract_character(pattern, i)
    print("\n".join(char_block))
    print()

for i in range(0, 35, 7):
    char_block = extract_character(pattern, i, start_row=8)
    print("\n".join(char_block))
    print()