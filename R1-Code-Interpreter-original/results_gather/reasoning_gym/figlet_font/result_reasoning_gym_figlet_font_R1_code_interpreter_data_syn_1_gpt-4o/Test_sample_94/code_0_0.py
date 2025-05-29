from pyfiglet import Figlet

# Create a Figlet object with the desired font style
figlet = Figlet(font='standard')

# Generate the figlet font for the word "figlet"
output = figlet.renderText('figlet')

# Print the output
print(output)