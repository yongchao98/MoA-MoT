def find_hook_flash_champions():
    """
    Analyzes which champions can perform a hook-flash to extend range.

    A hook-flash is defined as casting a hook and then flashing backward
    during the cast animation to increase the hook's effective range.
    """
    champions = {
        "Blitzcrank": False,
        "Nautilus": False,
        "Pyke": False,
        "Thresh": True  # Thresh can flash during his Q wind-up to change its origin point.
    }

    capable_champions = [name for name, can_do it in champions.items() if can_do]

    if not capable_champions:
        result = "None"
    else:
        result = ",".join(capable_champions)

    print(result)

find_hook_flash_champions()