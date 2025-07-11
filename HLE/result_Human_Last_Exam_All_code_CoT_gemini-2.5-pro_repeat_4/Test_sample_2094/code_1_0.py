def find_hook_flash_champions():
    """
    This function identifies and prints the champions from a predefined list
    who can perform a 'hook-flash' maneuver as described.
    """
    # List of champions to check
    all_champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # Champions who can perform the hook-flash
    # Blitzcrank, Nautilus, and Pyke can flash during their Q cast/charge time.
    # Thresh cannot flash during his Q wind-up.
    capable_champions = ["Blitzcrank", "Nautilus", "Pyke"]

    # Print the result as a comma-separated string
    result = ",".join(capable_champions)
    print(result)

find_hook_flash_champions()