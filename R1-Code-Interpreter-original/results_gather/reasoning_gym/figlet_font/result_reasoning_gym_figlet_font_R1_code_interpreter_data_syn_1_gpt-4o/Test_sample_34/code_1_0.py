import pyfiglet

# Define a simpler text to be converted to figlet font
text = "HELLO WORLD"

# Generate the figlet font
figlet_font = pyfiglet.figlet_format(text)

# Print the result
print(figlet_font)