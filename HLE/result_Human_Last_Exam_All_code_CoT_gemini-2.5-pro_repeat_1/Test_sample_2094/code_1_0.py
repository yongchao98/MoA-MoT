def find_hook_flash_champions():
    """
    Analyzes which champions can perform a 'hook-flash' as described:
    Hooking in one direction and flashing in the opposite direction to increase range.
    """
    champions = {
        "Blitzcrank": False,
        "Nautilus": False,
        "Pyke": False,
        "Thresh": True  # Thresh can flash during his Q wind-up, effectively extending its range from his initial cast position.
    }

    capable_champions = [name for name, can_do it in champions.items() if can_do]

    if not capable_champions:
        print("None")
    else:
        print(",".join(capable_champions))

find_hook_flash_champions()