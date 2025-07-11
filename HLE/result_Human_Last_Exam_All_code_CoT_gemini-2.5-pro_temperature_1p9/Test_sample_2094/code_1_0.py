def solve():
    """
    This function identifies which champions from a given list can perform a hook-flash combo
    to extend the hook's range.

    The champions to check are: Blitzcrank, Nautilus, Pyke, Thresh.

    - Thresh can cast his hook (Q) and flash during the wind-up animation. The hook will
      then originate from his new position, effectively extending its range.
    - Blitzcrank, Nautilus, and Pyke cannot use flash after casting their hook to extend
      its range. Their hooks originate from the location where the cast was finalized.

    Therefore, only Thresh fits the description.
    """
    answer = "Thresh"
    print(answer)

solve()
<<<Thresh>>>