# Define the ASCII art for each letter in the "standard" style
ascii_art_letters = {
    'H': [
        " _    _ ",
        "| |  | |",
        "| |__| |",
        "|  __  |",
        "| |  | |",
        "|_|  |_|"
    ],
    'E': [
        " ______ ",
        "|  ____|",
        "| |__   ",
        "|  __|  ",
        "| |____ ",
        "|______|"
    ],
    'L': [
        " _      ",
        "| |     ",
        "| |     ",
        "| |     ",
        "| |____ ",
        "|______|"
    ],
    'O': [
        "  ____  ",
        " / __ \ ",
        "| |  | |",
        "| |  | |",
        "| |__| |",
        " \____/ "
    ]
}

# The input ASCII art
input_art = [
    ".___  ____   ____ __________  ___________   _________   _________ ___________",
    "|   | \   \ /   / \______   \ \_   _____/  /   _____/  /   _____/ \_   _____/",
    "|   |  \   Y   /   |       _/  |    __)_   \_____  \   \_____  \   |    __)_ ",
    "|   |   \     /    |    |   \  |        \  /        \  /        \  |        \\",
    "|___|    \___/     |____|_  / /_______  / /_______  / /_______  / /_______  /",
    "                          \/          \/          \/          \/          \/ "
]

# Function to match input art to known letters
def match_ascii_art(input_art, ascii_art_letters):
    # Split the input art into columns for each letter
    columns = [input_art[i][j:j+9] for j in range(0, len(input_art[0]), 9) for i in range(6)]
    # Match each column to a letter
    result = []
    for col in columns:
        for letter, art in ascii_art_letters.items():
            if all(col[i] == art[i] for i in range(6)):
                result.append(letter)
                break
    return ''.join(result)

# Decode the input art
decoded_word = match_ascii_art(input_art, ascii_art_letters)
print(decoded_word)