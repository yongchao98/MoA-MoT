def find_hook_flash_champions():
    """
    This function determines which champions from a given list can perform a
    specific range-extending hook-flash maneuver.

    The maneuver is defined as casting a hook in one direction and flashing
    in the opposite direction during the cast animation to increase the hook's total range.

    - Blitzcrank, Nautilus, Pyke: Flashing repositions the origin of the hook,
      but flashing backward while hooking forward does not extend the range.
    - Thresh: Can aim his hook backward, flash forward during the wind-up, and the hook
      will travel backward from his new location, effectively extending its range.
      This matches the description.
    """
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    capable_champions = ["Thresh"]
    
    # The final answer is the list of champions who can perform the described mechanic.
    final_answer = ", ".join(capable_champions)
    
    print(final_answer)

find_hook_flash_champions()