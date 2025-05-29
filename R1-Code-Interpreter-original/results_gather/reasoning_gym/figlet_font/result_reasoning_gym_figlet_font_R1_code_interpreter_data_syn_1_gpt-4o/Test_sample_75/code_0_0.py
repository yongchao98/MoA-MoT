# Define the figlet text as a list of strings, each representing a line
figlet_lines = [
    " __   .__.  .  .   __.  .___.  .__.  .___.",
    "/  `  |  |  |\\ |  (__     |    [__]    |  ",
    "\\__.  |__|  | \\|  .__)    |    |  |    |  "
]

# Define the expected word
expected_word = "ASK"

# Print the figlet lines
for line in figlet_lines:
    print(line)

# Print the expected word
print("\nThe word represented by the figlet font is:", expected_word)