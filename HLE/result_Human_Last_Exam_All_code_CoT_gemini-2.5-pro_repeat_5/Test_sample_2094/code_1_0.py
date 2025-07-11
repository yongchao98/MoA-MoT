def solve_hook_flash_champions():
    """
    This function determines which of the listed champions can perform a "hook-flash"
    by buffering their hook ability with the Flash summoner spell.

    A hook-flash works if an ability has a cast time during which Flash can be used,
    causing the ability to originate from the new location.

    - Blitzcrank (Q): Has a cast time. Can Q-Flash.
    - Nautilus (Q): Has a cast time. Can Q-Flash.
    - Pyke (Q): Can charge Q, Flash, then release. Can Q-Flash.
    - Thresh (Q): Has a wind-up animation. Can Q-Flash.

    Therefore, all listed champions can perform this action.
    """
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    result = ",".join(champions)
    print(result)

solve_hook_flash_champions()