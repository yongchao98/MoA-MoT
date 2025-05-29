# Define the ASCII art as a list of strings
ascii_art = [
    "##                                   ##                     ",
    "##                                   ##                     ",
    "##                                   ##                     ",
    "#####                      ###    ########  ###    ######## ",
    "#####                      ####   ########  ####   ######## ",
    "##                           ###              ###           ",
    "##     ########  ########     ##               ##           ",
    "##     ########  ########     ##               ##           "
]

# Function to interpret the ASCII art
def interpret_ascii_art(ascii_art):
    # Define the patterns for letters
    patterns = {
        "H": [
            "##   ##",
            "##   ##",
            "##   ##",
            "######",
            "######",
            "##   ##",
            "##   ##",
            "##   ##"
        ],
        "I": [
            "######",
            "######",
            "  ##  ",
            "  ##  ",
            "  ##  ",
            "  ##  ",
            "######",
            "######"
        ]
    }
    
    # Check each pattern against the ASCII art
    for letter, pattern in patterns.items():
        if all(line.strip() in pattern for line in ascii_art):
            return letter
    return "Unknown"

# Interpret the ASCII art
word = interpret_ascii_art(ascii_art)
print(word)