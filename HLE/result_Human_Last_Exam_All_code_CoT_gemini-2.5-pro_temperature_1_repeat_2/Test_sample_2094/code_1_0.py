def find_hook_flash_champions():
    """
    This function identifies which champions from a given list can perform a
    hook-flash combo to extend the hook's range.

    The champions to check are: Blitzcrank, Nautilus, Pyke, Thresh.

    In League of Legends, the mechanic of casting a hook and then Flashing
    to have the hook originate from the new location is specific.

    - Blitzcrank's hook (Q) origin is set when cast.
    - Nautilus's hook (Q) origin is set when cast.
    - Pyke's hook (Q) origin is set when released.
    - Thresh's hook (Q) has a wind-up, during which he can Flash,
      and the hook will originate from his new position. This allows him
      to flash away from the target to extend the hook's effective range.

    The function will filter the list and print the name of the champion(s)
    who can perform this action.
    """
    champions_in_question = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # This list contains the champions known to be able to perform the described mechanic.
    capable_champions_list = ["Thresh"]

    # Find the intersection of the two lists.
    result = [champ for champ in champions_in_question if champ in capable_champions_list]

    if not result:
        print("None")
    else:
        # The join function handles the case of one or more champions.
        print(",".join(result))

find_hook_flash_champions()
<<<Thresh>>>