def find_hook_flash_champions():
    """
    This function analyzes a list of champions to determine which can perform
    a specific "hook-flash" combo to extend the ability's range.

    The analysis is based on established game mechanics:
    - Blitzcrank, Nautilus, and Pyke cannot extend their hook's range by flashing
      in the opposite direction after the cast has started. Their hooks are standard
      projectiles whose paths are locked.
    - Thresh is the only champion on the list who can flash during his Q's
      cast animation to change the hook's origin point relative to its trajectory,
      thereby extending its effective range.
    """
    
    # The list of champions capable of this specific mechanic.
    capable_champions = ["Thresh"]

    # If the list is empty, the result should be 'None'.
    # Otherwise, join the names with a comma.
    if not capable_champions:
        result = "None"
    else:
        result = ",".join(capable_champions)
        
    print(result)

find_hook_flash_champions()