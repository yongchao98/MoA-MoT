def solve_hook_flash_query():
    """
    This function determines which champions from a predefined list can perform
    a hook-flash combo to extend the ability's range.

    The analysis is based on in-game mechanics for each champion:
    - Blitzcrank: Cannot. His hook's origin is fixed immediately upon casting.
    - Nautilus: Cannot. His hook's origin is also fixed immediately upon casting.
    - Pyke: Can. He can charge his hook, Flash, and then release it from his new location.
    - Thresh: Can. He can Flash during his hook's wind-up animation, causing the hook to originate from his new location.
    """
    
    # List of champions who can perform the specified hook-flash combo.
    capable_champions = ["Pyke", "Thresh"]

    # Join the list elements into a single comma-separated string for the final answer.
    # The problem asks to output each name in the final result.
    # We will print the final string: "Pyke,Thresh"
    final_answer = ",".join(capable_champions)
    
    print(final_answer)

solve_hook_flash_query()