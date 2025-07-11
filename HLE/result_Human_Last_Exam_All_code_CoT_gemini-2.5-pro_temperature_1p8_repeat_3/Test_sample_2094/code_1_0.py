def find_hook_flash_champions():
    """
    Identifies which champions can perform a hook-flash by flashing
    in the opposite direction of the hook to extend its range.
    """
    # Based on in-game mechanics as of Season 14 (2024), only Thresh can
    # exploit a specific timing window to extend his hook range by flashing backwards.
    # The other champions' hooks originate from their post-flash position, which
    # would shorten the hook's effective range if they flashed backwards.
    champion_mechanics = {
        'Blitzcrank': False,
        'Nautilus': False,
        'Pyke': False,
        'Thresh': True
    }

    capable_champions = [name for name, can_perform_mechanic in champion_mechanics.items() if can_perform_mechanic]

    if not capable_champions:
        print("None")
    else:
        print(",".join(capable_champions))

find_hook_flash_champions()