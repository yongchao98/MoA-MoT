# Define the figlet font text
figlet_text = """
 #####      ###   ##  ##       ##    ####      ####       ###   
      #   ##  ##  ##  ##       ##   #    #    #    #    ##   #  
 ##  ##  #    ##  ##   ##     ##     # ##      # ##    #        
 ## ##    ######  ##   ##     ##     ###       ###      ###     
### #     #   #   ##  ##      ##   ##  ##    ##  ##     #       
 ##       #  ##    ####       ##        ##        ##    #   #   
 ##       #  ##    ###        ##    #######   #######   ####    
"""

# Define a mapping of figlet font patterns to standard text characters
figlet_to_text = {
    "#####": "H",
    "###": "E",
    "##  ##": "L",
    "####": "O"
}

# Decode the figlet font text
decoded_text = ""
for line in figlet_text.splitlines():
    for pattern, char in figlet_to_text.items():
        if pattern in line:
            decoded_text += char
            break

# Print the decoded text
print(decoded_text)