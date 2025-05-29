import pyfiglet

# Create a figlet object with a specific font
figlet = pyfiglet.Figlet(font='block')

# Generate ASCII art for the text "DEI GNN"
ascii_art = figlet.renderText("DEI GNN")

# Print the result
print(ascii_art)