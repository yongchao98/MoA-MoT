# Define the ASCII art as a list of strings
ascii_art = [
    "       ",
    "       ",
    "#######  ###      #####   ######   ",
    "###      ###     ###  ##  ###  ##  ",
    "######   ###     ###  ##  ###  ##  ",
    "###      ###     ###  ##  ###  ##  ",
    "###      ###     ###  ##  ###  ##  ",
    "#######  ######   #####   ######   "
]

# Function to print the ASCII art
def print_ascii_art(art):
    for line in art:
        print(line)

# Print the ASCII art
print_ascii_art(ascii_art)