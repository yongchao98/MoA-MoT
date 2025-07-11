def solve_hook_flash_champions():
    """
    Identifies which of the specified champions can perform a hook-flash.

    A hook-flash is a mechanic where a champion uses the Flash summoner spell
    during the cast time of their hook ability. This causes the hook to originate
    from the new, post-Flash location, effectively extending its threat range.
    """
    # The list of champions to check, as provided in the question.
    champions_in_question = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # My knowledge base: In League of Legends, all the listed champions have a
    # cast time on their hook ability, which allows Flash to be used during it.
    # Therefore, all of them can perform the hook-flash.
    capable_champions = []
    for champ in champions_in_question:
        # All listed champions can perform this mechanic.
        capable_champions.append(champ)

    # Format the result for printing.
    if not capable_champions:
        result = "None"
    else:
        # Join the list of names with a comma and space.
        result = ", ".join(capable_champions)

    print(result)

solve_hook_flash_champions()
<<<Blitzcrank, Nautilus, Pyke, Thresh>>>