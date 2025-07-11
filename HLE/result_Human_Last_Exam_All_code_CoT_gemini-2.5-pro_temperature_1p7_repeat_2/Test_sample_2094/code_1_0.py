def solve_hook_flash_riddle():
    """
    This function determines which champions from a given list can perform a "hook-flash"
    to extend their hook's range by flashing in the opposite direction of the cast.
    """
    
    # The list of champions to check.
    champions_to_check = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # A dictionary representing the knowledge about each champion's mechanics.
    # True if the champion can buffer their hook cast with Flash to extend its range.
    # False otherwise.
    mechanic_ability = {
        "Blitzcrank": False,
        "Nautilus": False,
        "Pyke": True,
        "Thresh": True
    }
    
    # List to store the names of champions who can perform the maneuver.
    capable_champions = []
    
    # Iterate through the list of champions and check their capability.
    for champion in champions_to_check:
        if mechanic_ability.get(champion, False):
            capable_champions.append(champion)

    # Format the final output.
    if not capable_champions:
        print("None")
    else:
        # Join the list of capable champions with a comma.
        result = ",".join(capable_champions)
        print(result)

solve_hook_flash_riddle()