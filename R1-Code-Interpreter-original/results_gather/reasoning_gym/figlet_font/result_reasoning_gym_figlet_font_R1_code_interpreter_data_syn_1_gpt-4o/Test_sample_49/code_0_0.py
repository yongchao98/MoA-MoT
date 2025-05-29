import pyfiglet

# Define the text to be converted into figlet font
text = "The User asks a question, and you solve it."

# Generate the figlet font
figlet_font = pyfiglet.figlet_format(text)

# Print the result
print(figlet_font)