def solve_champion_hook_flash():
    """
    This function identifies which of the listed champions can perform a 'hook-flash'
    to extend their hook's range.
    """
    champions_list = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # All four champions can perform a version of this mechanic.
    capable_champions = champions_list
    
    # Print the names separated by a comma.
    # The output shows each champion that makes up the final answer.
    print(",".join(capable_champions))

solve_champion_hook_flash()