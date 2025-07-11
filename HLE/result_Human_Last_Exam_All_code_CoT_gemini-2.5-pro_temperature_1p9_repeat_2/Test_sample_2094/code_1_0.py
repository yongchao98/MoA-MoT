def find_hook_flash_champions():
    """
    Analyzes which champions can perform a hook-flash by flashing
    in the opposite direction to extend the hook's range.
    """
    champions_data = {
        'Blitzcrank': {
            'can_extend_hook_by_reverse_flash': False,
            'notes': "Flashing during Q cast time changes the hook's origin. Flashing backward shortens its effective range."
        },
        'Nautilus': {
            'can_extend_hook_by_reverse_flash': False,
            'notes': "Similar to Blitzcrank, his hook originates from his post-flash position."
        },
        'Pyke': {
            'can_extend_hook_by_reverse_flash': False,
            'notes': "His hook is released from his post-flash position. Flashing backward shortens the throw."
        },
        'Thresh': {
            'can_extend_hook_by_reverse_flash': True,
            'notes': "Due to a specific interaction with his Q's wind-up, flashing backward at the right moment can extend the hook's effective range beyond its normal cast indicator."
        }
    }

    capable_champions = []
    for champion, data in champions_data.items():
        if data['can_extend_hook_by_reverse_flash']:
            capable_champions.append(champion)

    if not capable_champions:
        print("None")
    else:
        print(",".join(capable_champions))

find_hook_flash_champions()