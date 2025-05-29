# Define the figlet font text as a list of strings
figlet_font = [
    "  oo_           .-.           c  c     \\\    ///               \\\  /// ",
    " /  _)-<      c(O_O)c         (OO)     ((O)  (O))       /)      ((O)(O)) ",
    " \\__ `.      ,'.---.`,      ,'.--.)     | \\  / |      (o)(O)     | \\ ||  ",
    "    `. |    / /|_|_|\\ \\    / //_|_\\     ||\\\\//||       //\\\\      ||\\\\||  ",
    "    _| |    | \\_____/ |    | \\___       || \\/ ||      |(__)|     || \\ |  ",
    " ,-'   |    '. `---' .`    '.    )      ||    ||      /,-. |     ||  ||  ",
    "(_..--'       `-...-'        `-.'      (_/    \\_)    -'   ''    (_/  \\_) "
]

# Function to interpret the figlet font
def interpret_figlet(figlet_font):
    # Join the lines to form the complete message
    message = "\n".join(figlet_font)
    return message

# Print the interpreted message
print(interpret_figlet(figlet_font))