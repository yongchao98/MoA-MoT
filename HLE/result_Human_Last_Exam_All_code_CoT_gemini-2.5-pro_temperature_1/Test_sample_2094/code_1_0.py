def find_hook_flash_champs():
    """
    Identifies which champions from a predefined list can perform a hook-flash.
    The list includes Blitzcrank, Nautilus, Pyke, and Thresh.
    All of them are capable of this mechanic.
    """
    champions = [
        "Blitzcrank",
        "Nautilus",
        "Pyke",
        "Thresh"
    ]
    
    # Joining the list of champions with a comma and space
    result = ", ".join(champions)
    
    # Printing the final result
    print(result)

find_hook_flash_champs()