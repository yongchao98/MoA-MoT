import pyfiglet

# Define the text to be converted to figlet font
text = " ____      __      _  _      __   \n(  _ \\    /__\\    ( \\( )    /__\\  \n ) _ <   /(__)\\    )  (    /(__)\\ \n(____/  (__)(__)  (_)_\\_)  (__)(__)\n"

# Use pyfiglet to convert the text
figlet_font = pyfiglet.figlet_format(text)

# Print the result
print(figlet_font)