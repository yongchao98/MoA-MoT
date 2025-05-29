# Define the pattern for each letter in a 7x8 grid
letter_patterns = {
    'H': [
        "## ### #",
        "## # # #",
        "## # # #",
        "###### #",
        "## # # #",
        "## # # #",
        "## ### #"
    ],
    'E': [
        "###### #",
        "## ### #",
        "## # # #",
        "###### #",
        "## # # #",
        "## ### #",
        "###### #"
    ],
    'L': [
        "## ### #",
        "## # # #",
        "## # # #",
        "## # # #",
        "## # # #",
        "## ### #",
        "###### #"
    ],
    'O': [
        "###### #",
        "## ### #",
        "## # # #",
        "## # # #",
        "## # # #",
        "## ### #",
        "###### #"
    ]
}

# Input pattern
input_pattern = [
    "###### #  ###### #  ###### #  ###### #  ###### #  ## ### #  ###### # ",
    "## ### #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## ### # ",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## # # # ",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ###### # ",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #   # ### # ",
    "## ### #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## ### # ",
    "###### #  ###### #   ### # #  ###### #  ## ### #  ###### #  ###### # ",
    " # # # #   # # # #   # # # #   # # # #   # # # #   # # # #   # # # # "
]

# Split the input into blocks of 7 lines each
letter_blocks = []
for i in range(0, len(input_pattern[0]), 9):  # 9 includes the space between letters
    block = [line[i:i+8] for line in input_pattern[:7]]
    letter_blocks.append(block)

# Function to identify the letter from a block
def identify_letter(block):
    for letter, pattern in letter_patterns.items():
        if block == pattern:
            return letter
    return '?'

# Identify each letter
word = ''.join(identify_letter(block) for block in letter_blocks)

print(word)