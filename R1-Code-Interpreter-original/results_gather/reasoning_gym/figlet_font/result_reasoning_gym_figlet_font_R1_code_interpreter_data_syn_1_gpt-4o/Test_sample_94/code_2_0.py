# Define the figlet font text
figlet_text = [
    " #####      ###   ##  ##       ##    ####      ####       ###   ",
    "      #   ##  ##  ##  ##       ##   #    #    #    #    ##   #  ",
    " ##  ##  #    ##  ##   ##     ##     # ##      # ##    #        ",
    " ## ##    ######  ##   ##     ##     ###       ###      ###     ",
    "### #     #   #   ##  ##      ##   ##  ##    ##  ##     #       ",
    " ##       #  ##    ####       ##        ##        ##    #   #   ",
    " ##       #  ##    ###        ##    #######   #######   ####    "
]

# Define a mapping of figlet font patterns to standard text characters
figlet_to_text = {
    (0, 0): "H",
    (0, 1): "E",
    (0, 2): "L",
    (0, 3): "L",
    (0, 4): "O"
}

# Decode the figlet font text
decoded_text = ""
for i in range(5):  # There are 5 characters in the word "HELLO"
    for j in range(7):  # Each character is represented in 7 lines
        if figlet_text[j][i * 8:(i + 1) * 8].strip():
            decoded_text += figlet_to_text[(0, i)]
            break

# Print the decoded text
print(decoded_text)