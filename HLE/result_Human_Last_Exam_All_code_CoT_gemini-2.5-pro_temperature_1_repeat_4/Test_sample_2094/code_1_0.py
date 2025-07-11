def find_hook_flash_champions():
    """
    This function identifies which of the selected champions can perform a "hook-flash".

    A "hook-flash" is a mechanic where a champion uses the Flash summoner spell
    during the cast time or wind-up of their hook ability. This causes the hook
    to originate from their new location after the flash, effectively extending
    the hook's range or changing its angle to surprise an opponent.

    Based on the game mechanics in Season 14 (2024):
    - Blitzcrank (Q - Rocket Grab): Can flash during the cast time.
    - Nautilus (Q - Dredge Line): Can flash during the cast time.
    - Pyke (Charged Q - Bone Skewer): Can flash while charging his hook before release.
    - Thresh (Q - Death Sentence): Can flash during the wind-up animation.

    All four champions can perform this technique.
    """
    capable_champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    print(",".join(capable_champions))

find_hook_flash_champions()