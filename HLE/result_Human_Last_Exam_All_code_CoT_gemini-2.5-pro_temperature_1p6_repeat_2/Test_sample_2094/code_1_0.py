def solve_hook_flash_question():
    """
    Analyzes which champions can perform a hook-flash as described:
    hooking forward and flashing backward to increase the hook's range.
    """
    champions = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]
    capable_champions = []

    # In League of Legends, ability buffering with Flash allows some champions
    # to cast an ability and then Flash during the animation to alter its effect.
    #
    # - Blitzcrank/Nautilus: Their hooks are projectiles cast from their location.
    #   Flashing doesn't extend the range in the way described.
    # - Pyke: He can Flash during his hook's channel, but the hook is then thrown
    #   from his new location, not extending the range if he Flashes backward.
    # - Thresh: He can cast his hook (Q), and during the wind-up animation, Flash
    #   backward. The hook still travels in the original direction but from his new,
    #   further-back position, effectively increasing its maximum range.
    #
    # Therefore, only Thresh fits the specific description.
    
    for champion in champions:
        if champion == "Thresh":
            capable_champions.append(champion)

    if not capable_champions:
        result = "None"
    else:
        result = ",".join(capable_champions)

    print("The champion(s) that can perform the described hook-flash are:")
    print(result)
    print(f"<<<{result}>>>")

solve_hook_flash_question()