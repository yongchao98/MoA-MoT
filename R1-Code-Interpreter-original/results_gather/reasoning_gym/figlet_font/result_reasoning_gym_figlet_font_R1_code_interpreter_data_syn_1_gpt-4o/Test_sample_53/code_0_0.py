# Define the ASCII art pattern
ascii_art = [
    "   ##      ###     ####       ##      ###    ####      ###   ",
    " #   #   ##  ##   #    #    #   #   ##  ##  #    #   ##  ##  ",
    "##      #    ##    # ##    ##      #    ##  #   #   #    ##  ",
    "##       ######    ###     ##       ######  #####    ######  ",
    "##   #   #   #   ##  ##    ##   #   #   #   #   ##   #   #   ",
    "##  #    #  ##        ##   ##  #    #  ##   #   #    #  ##   ",
    " ###     #  ##    #######   ###     #  ##   #  #     #  ##   ",
    "                                               #             "
]

# Define a function to interpret the ASCII art
def interpret_ascii_art(ascii_art):
    # Define a dictionary mapping ASCII patterns to letters
    ascii_to_letter = {
        " ###  #   #  ###  ###   #   #  #####  ###   #####  ###  ": "A",
        "####  #   #  ####  ###  #   #  #####  ###   #####  ###  ": "B",
        " ###  #   #  #     #    #   #  #      #     #     #     ": "C",
        "####  #   #  ####  ###  #   #  #####  ###   #####  ###  ": "D",
        "##### #     ##### #     #     #####  #     #####  #     ": "E",
        "##### #     ##### #     #     #####  #     #     #      ": "F",
        " ###  #     #     #  ## #   #  ###   #   #  ###   ###   ": "G",
        "#   # #   # ##### #   # #   # #####  #   # #####  #   # ": "H",
        " ###   #     #     #     #     #     #     #     ###   ": "I",
        "  ###   #     #     #     #     #     #   #  ###   ###  ": "J",
        "#   # #  #  ###   #  #  #   # #   #  #   # #   #  #   # ": "K",
        "#     #     #     #     #     #     #     #####  #####  ": "L",
        "###   # #   # #   # #   # #   # #   # #   # #   # #   # ": "M",
        "#   # ##  # #  ## #   # #   # #   # #   # #   # #   #  ": "N",
        " ###  #   # #   # #   # #   # #   # #   #  ###   ###   ": "O",
        "####  #   # #   # ####  #     #     #     #     #      ": "P",
        " ###  #   # #   # #   # #   # #  ##  ###   ###   ###   ": "Q",
        "####  #   # #   # ####  #  #  #   # #   # #   #  #   # ": "R",
        " ###  #   # #     ###    ###     #  #   #  ###   ###   ": "S",
        "#####  #     #     #     #     #     #     #     #     ": "T",
        "#   # #   # #   # #   # #   # #   # #   #  ###   ###   ": "U",
        "#   # #   # #   # #   # #   #  # #    # #    #     #   ": "V",
        "#   # #   # #   # #   # # # # # # # # # # # # # # # #  ": "W",
        "#   #  # #   # #    #     #     #    # #   # #   #   # ": "X",
        "#   #  # #   # #    #     #     #     #     #     #    ": "Y",
        "#####    #    #    #    #    #    #    #    ##### ##### ": "Z"
    }
    
    # Extract the letters from the ASCII art
    letters = []
    for i in range(0, len(ascii_art[0]), 7):
        letter_pattern = "".join([line[i:i+5] for line in ascii_art])
        if letter_pattern in ascii_to_letter:
            letters.append(ascii_to_letter[letter_pattern])
        else:
            letters.append("?")
    
    return "".join(letters)

# Interpret the ASCII art
word = interpret_ascii_art(ascii_art)
print(word)