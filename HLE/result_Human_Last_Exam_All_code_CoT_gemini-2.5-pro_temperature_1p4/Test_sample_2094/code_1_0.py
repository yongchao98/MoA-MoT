def find_hook_flash_champions():
    """
    Analyzes which of the given champions can perform a specific hook-flash mechanic.

    The mechanic is defined as hooking in one direction and flashing in the opposite
    direction to effectively increase the hook's range relative to the champion's
    final position.
    """
    
    # List of champions to evaluate
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # In League of Legends, an ability's origin point is key.
    # For Blitzcrank, Nautilus, and Pyke, flashing during the cast/charge time
    # moves the ability's origin point to their new location. Flashing backwards
    # would therefore DECREASE the hook's range toward a target.
    #
    # Thresh is the unique case. The origin of his Q (Death Sentence) is locked
    # at the START of his cast animation. He can then flash backwards during the
    # wind-up, and the hook will still fire from the original location. This
    # matches the description of increasing the hook's range.
    
    can_perform_mechanic = ["Thresh"]
    
    # Filter the original list to find the correct champions
    result = [champion for champion in champions if champion in can_perform_mechanic]
    
    # Format the final answer
    if not result:
        final_answer = "None"
    else:
        final_answer = ",".join(result)
        
    print(f"Champions who can perform the specified hook-flash: {final_answer}")
    print(f"<<<{final_answer}>>>")

find_hook_flash_champions()