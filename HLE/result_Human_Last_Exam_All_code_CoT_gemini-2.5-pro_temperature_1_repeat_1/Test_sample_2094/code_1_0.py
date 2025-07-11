def find_hook_flash_champions():
    """
    Identifies which champions from a predefined list can perform a hook-flash.

    A hook-flash is defined as casting a hook and then flashing in the opposite
    direction during the cast time to extend the ability's effective range.
    Based on League of Legends mechanics as of Season 14 (2024), all listed
    champions (Blitzcrank, Nautilus, Pyke, Thresh) can perform this action.
    """
    champions = [
        "Blitzcrank",
        "Nautilus",
        "Pyke",
        "Thresh"
    ]
    
    # All champions in the list can perform this mechanic.
    # We will join their names with a comma for the final output.
    result = ",".join(champions)
    
    print(result)

find_hook_flash_champions()