def solve_hook_flash_question():
    """
    Analyzes which of the given champions can perform a "hook-flash"
    by hooking in one direction and flashing in the opposite direction
    to increase the hook's effective range.

    This analysis is based on in-game mechanics for League of Legends as of Season 14 (2024).

    - Blitzcrank: Cannot. Flash changes the hook's origin point.
    - Nautilus: Cannot. Flash changes the hook's origin point.
    - Pyke: Cannot. Flash changes the hook's origin point.
    - Thresh: Can. A specific timing of casting Q and flashing backwards makes the hook projectile travel farther.
    """
    
    champions_who_can = ["Thresh"]
    
    if not champions_who_can:
        result = "None"
    else:
        result = ",".join(champions_who_can)
        
    print(result)

solve_hook_flash_question()