def find_hook_flash_champs():
    """
    Identifies which champions from a given list can perform a hook-flash combo.

    The "hook-flash" is a mechanic where a champion uses the Flash summoner spell
    during the cast or channel time of their hook ability to change its origin point,
    effectively extending its range.

    The champions evaluated are: Blitzcrank, Nautilus, Pyke, and Thresh.
    """
    
    # List of champions who can perform a hook-flash to extend range.
    # Blitzcrank (Q - Rocket Grab): Can flash during the cast time.
    # Nautilus (Q - Dredge Line): Can flash during the cast time.
    # Pyke (Q - Bone Skewer): Can flash while charging the hook.
    # Thresh (Q - Death Sentence): Can flash during the wind-up animation.
    
    capable_champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # Print the names separated by a comma.
    print(",".join(capable_champions))

find_hook_flash_champs()