def solve_hook_flash_champions():
    """
    This function determines which champions from a given list can perform a hook-flash
    to extend the ability's range.

    The mechanic works for champions whose hook's origin point is determined
    at the end of the cast animation, allowing Flash to be used during the cast
    to reposition the origin.

    - Blitzcrank (Q - Rocket Grab): No. The hook's origin is set at the start of the cast.
    - Nautilus (Q - Dredge Line): No. The hook's origin is set at the start of the cast.
    - Pyke (Q - Bone Skewer): Yes. Can Flash during the charge-up, and the hook releases from the new position.
    - Thresh (Q - Death Sentence): Yes. Can Flash during the wind-up animation, and the hook releases from the new position.
    """

    champion_abilities = {
        'Blitzcrank': False,
        'Nautilus': False,
        'Pyke': True,
        'Thresh': True
    }

    capable_champions = []
    for champion, can_perform in champion_abilities.items():
        if can_perform:
            capable_champions.append(champion)

    if capable_champions:
        result = ",".join(capable_champions)
    else:
        result = "None"
        
    print(result)

solve_hook_flash_champions()
<<<Pyke,Thresh>>>