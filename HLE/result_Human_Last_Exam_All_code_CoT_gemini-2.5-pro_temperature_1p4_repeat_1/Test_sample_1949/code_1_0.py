import operator

def find_longest_rivers_in_aeneid():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    The script contains a dictionary of major rivers mentioned in the text
    and their approximate modern lengths. It sorts them to find the top three.
    """

    # A dictionary of prominent rivers mentioned in the Aeneid and their lengths in km.
    # Sources include mentions on Aeneas's shield (Book VIII), prophecies, and settings.
    rivers = {
        'Nile': 6650,      # Mentioned in Book VIII, associated with Egypt and Cleopatra.
        'Indus': 3180,     # Mentioned in Book VIII, representing the far East.
        'Danube': 2850,    # Mentioned as 'Hister' in Book VIII, on the shield.
        'Euphrates': 2800, # Mentioned in Book VIII on the shield, representing Parthia.
        'Ganges': 2525,    # Mentioned in Book VIII, symbolizing distant conquests.
        'Rhine': 1230,     # Mentioned in Book VIII, representing the Gauls.
        'Po': 652,         # Mentioned as 'Eridanus' in the underworld in Book VI.
        'Tiber': 406,      # The river of Rome, central to the latter half of the Aeneid.
        'Araxes': 1072,    # Mentioned on the shield in Book VIII.
        'Hebrus': 480       # A Thracian river mentioned in Book I.
    }

    # Sort the rivers by length in descending order.
    # The 'sorted' function returns a list of tuples (key, value).
    sorted_rivers = sorted(rivers.items(), key=operator.itemgetter(1), reverse=True)

    # Get the top three longest rivers.
    top_three_rivers = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    for i, (river, length) in enumerate(top_three_rivers):
        # The prompt asks to "output each number in the final equation!".
        # We will format it as a clear statement showing the river and its length.
        print(f"{i+1}. {river} = {length} km")

if __name__ == '__main__':
    find_longest_rivers_in_aeneid()