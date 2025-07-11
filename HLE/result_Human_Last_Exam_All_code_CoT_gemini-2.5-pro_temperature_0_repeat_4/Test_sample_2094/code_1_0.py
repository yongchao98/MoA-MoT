def find_hook_flash_champions():
    """
    This function identifies which champions from a predefined list can perform a
    hook-flash by casting their hook and then flashing in the opposite direction
    to extend the hook's range.

    Based on League of Legends game mechanics as of Season 14 (2024):
    - Blitzcrank, Nautilus, and Pyke fire a projectile with a fixed range.
      Flashing backward after the cast does not extend the hook's range.
    - Thresh's hook (Death Sentence) is tethered to his position. Flashing
      backward while the hook is in flight will extend its maximum range.
    """
    
    # The list of champions provided by the user.
    champion_pool = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    
    # The list of champions who can perform the described mechanic.
    capable_champions = ["Thresh"]

    # Filter the pool to find the capable champions. This is more for show,
    # as the logic is based on external game knowledge.
    result = [champ for champ in champion_pool if champ in capable_champions]

    if not result:
        print("None")
    else:
        # Join the names with a comma for the final output.
        print(",".join(result))

find_hook_flash_champions()