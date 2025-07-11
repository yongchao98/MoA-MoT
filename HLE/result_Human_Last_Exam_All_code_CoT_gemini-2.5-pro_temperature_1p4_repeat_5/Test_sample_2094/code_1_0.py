def solve():
    """
    This function identifies which champions from a given list can perform a
    "hook-flash" combo, where they flash backward during the cast time to
    extend the hook's range.
    """

    # The list of champions provided by the user.
    champions_to_check = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # Champions who can perform the specified hook-flash mechanic.
    # For Blitzcrank, Nautilus, and Pyke, the hook's origin is updated to their
    # new position after flashing during the cast animation.
    # For Thresh, the hook's origin is locked at the start of the cast and is
    # not updated by Flash.
    capable_champions = ["Blitzcrank", "Nautilus", "Pyke"]

    # Print the names of the capable champions, separated by a comma.
    result = ",".join(capable_champions)
    print(result)

solve()
<<<Blitzcrank,Nautilus,Pyke>>>