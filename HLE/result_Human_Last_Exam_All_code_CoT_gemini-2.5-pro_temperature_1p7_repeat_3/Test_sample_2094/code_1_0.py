def solve_hook_flash_question():
    """
    This function determines which champions from a given list can perform a hook-flash,
    which is defined as casting the hook and then using Flash to change the hook's origin point
    during its cast/channel time to increase its effective range.

    - Thresh: Can Flash during his Q's wind-up animation. The hook originates from his new location.
    - Pyke: Can Flash while channeling his Q. The hook originates from his new location upon release.
    - Blitzcrank/Nautilus: Their hook origins are locked at the start of the cast, so Flash cannot be used to alter it mid-cast.
    """

    all_champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    capable_champions = []

    # Logic based on in-game mechanics for Season 14 (2024)
    for champion in all_champions:
        if champion == "Pyke" or champion == "Thresh":
            capable_champions.append(champion)

    if len(capable_champions) > 0:
        # Join the names with a comma as requested
        result = ",".join(capable_champions)
    else:
        result = "None"
    
    print(result)

solve_hook_flash_question()