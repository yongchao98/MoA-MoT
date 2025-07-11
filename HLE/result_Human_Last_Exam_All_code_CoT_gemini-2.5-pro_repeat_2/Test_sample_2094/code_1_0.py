def solve_champion_mechanic():
    """
    This function determines which champions from a given list can perform a specific
    "hook-flash" mechanic as described.
    """
    # The list of champions provided in the problem.
    champions_to_check = ["Blitzcrank", "Nautilus", "Pyke", "Thresh"]

    # A knowledge base mapping each champion to their ability to perform the hook-flash.
    # True means they can flash during the hook's cast animation to change its origin point.
    # False means they cannot.
    mechanic_ability = {
        "Blitzcrank": False,
        "Nautilus": False,
        "Pyke": False,
        "Thresh": True
    }

    # List to store the names of champions who can perform the mechanic.
    capable_champions = []

    # Iterate through the provided list and check our knowledge base.
    for champion in champions_to_check:
        if mechanic_ability.get(champion, False):
            capable_champions.append(champion)

    # Format the result based on whether any champions were found.
    if not capable_champions:
        final_answer = "None"
    else:
        final_answer = ", ".join(capable_champions)
    
    print(final_answer)

solve_champion_mechanic()
<<<Thresh>>>