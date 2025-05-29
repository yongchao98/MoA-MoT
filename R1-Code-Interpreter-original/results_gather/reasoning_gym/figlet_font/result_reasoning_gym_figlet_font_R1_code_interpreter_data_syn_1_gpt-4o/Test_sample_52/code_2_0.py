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
    # Adjust the width of each letter based on the ASCII art
    letter_widths = [7, 7, 7, 7]  # Assuming each letter is 7 characters wide
    start_indices = [0, 9, 18, 27, 36]  # Manually determined start indices for each letter
    for start_index in start_indices:
        letter = [line[start_index:start_index+7] for line in ascii_art]
        for char, pattern in letter_patterns.items():
            if letter == pattern:
                letters.append(char)
                break
    
    return ''.join(letters)

# Extract and print the word
word = extract_letters(ascii_art)
print(word)