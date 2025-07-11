def solve():
    """
    This function identifies and prints the names of the champions who can
    perform a hook-flash, which is defined as casting a hook and then flashing
    during the cast animation to change the hook's origin point and extend its range.
    """
    
    # List of champions to check
    champions_to_check = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # Champions who can perform the hook-flash mechanic as described
    capable_champions = ["Pyke", "Thresh"]
    
    # Print the result separated by a comma
    print(",".join(capable_champions))

solve()