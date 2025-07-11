def solve():
    """
    This function identifies League of Legends champions from a predefined list
    who can perform a "Dash-Flash" combo and prints their names.
    The Dash-Flash combo is defined as using Flash during a dash ability's
    animation to alter the ability's point of impact for its effect (damage/CC).
    """

    # List of champions who can perform the Dash-Flash combo as defined.
    # The combo involves buffering Flash during the dash to extend the ability's effect.
    champions_who_can_dash_flash = [
        "Gragas",      # E -> Flash
        "Vi",          # Q -> Flash
        "Jarvan",      # E -> Q -> Flash
        "Galio",       # E -> Flash
        "Orn",         # E -> Flash
        "Rakan",       # W -> Flash
        "Pyke",        # E -> Flash
        "Rell",        # W (dismount) -> Flash
        "Riven",       # Q3 -> Flash
        "Sejuani",     # Q -> Flash
        "Shen",        # E -> Flash
        "Tryndamere",  # E -> Flash
        "Urgot",       # E -> Flash
        "Yone",        # Q3 -> Flash
        "Zac",         # E -> Flash
        "Camille"      # E2 (Wall Dive) -> Flash
    ]

    # Join the list of names with ", " and print the result.
    result = ", ".join(champions_who_can_dash_flash)
    print(result)

solve()
# The final list is constructed based on the in-game mechanic where Flash can be
# used during the travel time of a dash ability to change the destination and
# apply the ability's crowd control or damage at the new location.
# Champions whose dashes cannot be buffered this way (e.g., Fizz, LeBlanc, Tristana) are excluded.
# <<<Gragas, Vi, Jarvan, Galio, Orn, Rakan, Pyke, Rell, Riven, Sejuani, Shen, Tryndamere, Urgot, Yone, Zac, Camille>>>