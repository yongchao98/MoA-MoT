def find_hook_flash_champions():
    """
    Analyzes a list of champions to see which can perform a reverse hook-flash
    to extend the ability's range.

    In League of Legends, this is a unique interaction primarily associated with Thresh.
    When Thresh casts his Q (Death Sentence), there's a wind-up animation. If he
    flashes backward during this wind-up, the hook's origin point moves back, but it
    still travels in the originally targeted direction. This effectively increases
    the hook's maximum range from where Thresh ends up.

    Other champions like Blitzcrank, Nautilus, and Pyke can Flash during their
    hook cast, but this is typically used to change the angle or surprise an
    opponent, not to extend the range by flashing backward.
    """
    
    # A dictionary representing the knowledge about each champion's mechanic.
    # True if they can perform the described reverse hook-flash for range.
    champion_mechanics = {
        'Blitzcrank': False,
        'Nautilus': False,
        'Pyke': False,
        'Thresh': True 
    }

    capable_champions = []
    # Iterate through the provided list of champions
    for champion in ['Blitzcrank', 'Nautilus', 'Pyke', 'Thresh']:
        if champion_mechanics.get(champion, False):
            capable_champions.append(champion)

    # Format the output according to the instructions
    if not capable_champions:
        print("None")
    else:
        # Join the list of capable champions with a comma
        result = ",".join(capable_champions)
        print(result)

find_hook_flash_champions()
<<<Thresh>>>