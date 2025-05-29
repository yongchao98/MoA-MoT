# Define the ASCII art as a list of strings
ascii_art = [
    "                              ##                       ",
    "                              ##                       ",
    " ##  ##    ####    #####     #####     ####     #####  ",
    " ##  ##   ##  ##   ##  ##     ##          ##   ##      ",
    " ##  ##   ######   ##  ##     ##       #####    #####  ",
    "  #####   ##       ##  ##     ## ##   ##  ##        ## ",
    "     ##    #####   ##  ##      ###     #####   ######  ",
    " #####                                                 "
]

# Define a function to interpret the ASCII art
def interpret_ascii_art(ascii_art):
    # Mapping of ASCII patterns to characters
    ascii_to_char = {
        " ##  ##    ####    #####     #####     ####     #####  ": "H",
        " ##  ##   ##  ##   ##  ##     ##          ##   ##      ": "E",
        " ##  ##   ######   ##  ##     ##       #####    #####  ": "L",
        "  #####   ##       ##  ##     ## ##   ##  ##        ## ": "L",
        "     ##    #####   ##  ##      ###     #####   ######  ": "O"
    }
    
    # Extract the characters from the ASCII art
    word = ""
    for line in ascii_art[2:7]:  # Only consider the lines that form the letters
        if line in ascii_to_char:
            word += ascii_to_char[line]
    
    return word

# Interpret the ASCII art
word = interpret_ascii_art(ascii_art)
print(word)