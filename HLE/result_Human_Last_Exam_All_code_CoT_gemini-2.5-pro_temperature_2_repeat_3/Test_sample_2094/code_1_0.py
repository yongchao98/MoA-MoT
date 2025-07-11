def find_hook_flash_champions():
    """
    Analyzes which champions can perform a specific "hook-flash" mechanic.
    
    The mechanic is defined as: hooking in one direction and flashing
    in the opposite direction to increase the hook's effective range.
    
    - Blitzcrank, Nautilus, Pyke: Flashing repositions the start of the hook.
      To extend range, they must flash *towards* the target. Flashing away shortens it.
      
    - Thresh: Can flash *after* the hook is in the air. Flashing away from the
      hook's direction increases the hook's maximum travel distance. This matches the criteria.
    """
    
    champions = {
        'Blitzcrank': False,
        'Nautilus': False,
        'Pyke': False,
        'Thresh': True 
    }
    
    capable_champions = [name for name, can_do_it in champions.items() if can_do_it]
    
    if not capable_champions:
        print("None")
    else:
        # The join function handles the comma separation as requested.
        print(",".join(capable_champions))

find_hook_flash_champions()