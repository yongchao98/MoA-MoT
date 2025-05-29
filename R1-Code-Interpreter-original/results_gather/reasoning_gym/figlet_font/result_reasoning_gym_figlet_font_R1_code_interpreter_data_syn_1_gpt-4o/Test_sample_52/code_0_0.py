# Define the ASCII art
ascii_art = [
    "######   #######    #  ##   #####   ######   #######  ##       ",
    "     ##            ###     ###  ##       ##                    ",
    "##   ##  ##       #######  ##   ##  ##   ##  ##       ##       ",
    "######   #####    ## # ##  ##   ##  ##   ##  #####    ##       ",
    "## ##    ##       ## # ##  ##   ##  ##  ##   ##       ##       ",
    "##  ##   ###  ##  ##   ##  ##  ###  ## ###   ###  ##  ##  ###  ",
    " ##  ##   #####    ##  ##   #####   ####      #####   ######   "
]

# Define a function to extract letters from the ASCII art
def extract_letters(ascii_art):
    # Known patterns for each letter
    letter_patterns = {
        'H': [
            "##   ##",
            "##   ##",
            "###### ",
            "##   ##",
            "##   ##",
            "##   ##",
            "##   ##"
        ],
        'E': [
            "###### ",
            "##     ",
            "##     ",
            "#####  ",
            "##     ",
            "##     ",
            "###### "
        ],
        'L': [
            "##     ",
            "##     ",
            "##     ",
            "##     ",
            "##     ",
            "##     ",
            "###### "
        ],
        'O': [
            " ##### ",
            "##   ##",
            "##   ##",
            "##   ##",
            "##   ##",
            "##   ##",
            " ##### "
        ]
    }
    
    # Extract letters from the ASCII art
    letters = []
    for i in range(0, len(ascii_art[0]), 8):  # Assuming each letter is 8 characters wide
        letter = [line[i:i+8] for line in ascii_art]
        for char, pattern in letter_patterns.items():
            if letter == pattern:
                letters.append(char)
                break
    
    return ''.join(letters)

# Extract and print the word
word = extract_letters(ascii_art)
print(word)