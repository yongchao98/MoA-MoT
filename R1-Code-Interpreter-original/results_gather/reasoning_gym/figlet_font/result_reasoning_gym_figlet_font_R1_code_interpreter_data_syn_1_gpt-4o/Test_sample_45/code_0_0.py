# Define the pattern as a list of strings
pattern = [
    "###### #  ###### #  ###### #  ###### #  ###### #  ## ### #  ###### #",
    "## ### #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## ### #",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## # # #",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ###### #",
    "## # # #  ## ### #   ### # #   ### # #  ## ### #  ## ### #   # ### #",
    "## ### #  ## ### #   ### # #   ### # #  ## ### #  ## ### #  ## ### #",
    "###### #  ###### #   ### # #  ###### #  ## ### #  ###### #  ###### #",
    " # # # #   # # # #   # # # #   # # # #   # # # #   # # # #   # # # #"
]

# Function to extract and print each letter
def extract_letters(pattern):
    # Number of lines per letter
    lines_per_letter = 8
    # Width of each letter (approximate)
    letter_width = 8
    
    # Number of letters
    num_letters = len(pattern[0]) // letter_width
    
    # Extract each letter
    for i in range(num_letters):
        print(f"Letter {i+1}:")
        for line in pattern:
            print(line[i*letter_width:(i+1)*letter_width])
        print("\n")

# Extract and print the letters
extract_letters(pattern)