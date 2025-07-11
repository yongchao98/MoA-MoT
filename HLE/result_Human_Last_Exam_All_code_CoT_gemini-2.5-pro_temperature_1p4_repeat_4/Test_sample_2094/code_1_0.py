def solve_hook_flash_query():
    """
    This function determines which champion(s) from a given list can perform
    a specific "hook-flash" mechanic.

    The mechanic is defined as: hooking a target, and then flashing in the
    opposite direction to increase the hook's effective range (or pull distance).

    Based on League of Legends mechanics (as of Season 14, 2024):
    - Blitzcrank: Cannot flash once the hook has landed.
    - Nautilus: Cannot flash once the hook has landed.
    - Pyke: Cannot flash once the hook has landed and the pull starts.
    - Thresh: CAN flash after the hook lands but before reactivation. Flashing
              away from the target slingshots them a greater distance. This is
              the specific interaction described.
    """
    champion_mechanics = {
        "Blitzcrank": False,
        "Nautilus": False,
        "Pyke": False,
        "Thresh": True
    }

    capable_champions = []
    for champion, can_perform_mechanic in champion_mechanics.items():
        if can_perform_mechanic:
            capable_champions.append(champion)

    if not capable_champions:
        final_answer = "None"
    else:
        final_answer = ",".join(capable_champions)

    print(final_answer)

solve_hook_flash_query()
<<<Thresh>>>