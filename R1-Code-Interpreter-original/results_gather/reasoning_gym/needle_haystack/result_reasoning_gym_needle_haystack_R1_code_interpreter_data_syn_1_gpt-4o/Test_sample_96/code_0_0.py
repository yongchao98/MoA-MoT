# Define the text
text = """
Marcin exults in tap dancing. Rohaan glories in acting. Lachlan finds fulfillment in digital art. 
Conal savors polishing the silver. Ralfs relishes roller coaster riding. Alistar pines for playing the oboe. 
Danyil treasures sandwich. Kames worships foxes. Sergio begrudges playing rugby. Henry loathes scrapbooking. 
Thumbiko revels in filmmaking. Stefano ignores hip-hop dancing. Jon-Paul pines for motorcycles. 
Azim exalts playing saxophone. Carlos is neutral toward rapping. Allan-Laiton appreciates kangaroos. 
Jeswin bemoans peacocks. Lennen eschews elegance. Christian adores playing badminton. 
Nasser is obsessed with buffaloes. Declyn lusts after crocheting. Makensie rejoices in playing tennis. 
Rui dotes urban exploration. Kristofer brushes off playing the ukulele. Pedram approves of vegetable soup. 
Kaydin basks in steak. Reis ridicules playing golf. Karim worships quantum physics. 
Saif is neutral toward chopping vegetables. Kiyonari prizes sailing. Mathias desires the color lavender. 
Justan overlooks playing field hockey. Hosea sneers at making coffee. Elyan adores swimming. 
Tane despises the color gray. Abaan hates machine learning. Corey-Jay scorns the color indigo. 
Tadhg puts up with the color black. Daksh execrates the color indigo. Azeem hates the color khaki.
"""

# Search for the person who relishes roller coaster riding
import re

# Use regex to find the name associated with "relishes roller coaster riding"
match = re.search(r'(\w+) relishes roller coaster riding', text)

# Print the name if found
if match:
    print(match.group(1))
else:
    print("No match found")