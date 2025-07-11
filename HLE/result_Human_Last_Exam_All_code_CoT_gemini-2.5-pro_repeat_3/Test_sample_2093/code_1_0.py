def find_dash_flash_champions():
    """
    Identifies champions from a predefined list who can perform a Dash-Flash combo.
    The combo is defined as using Flash during a dash ability's animation to
    change the location where the ability's primary effect (damage/CC) lands.
    """
    champions = [
        "Gragas",      # E -> Flash
        "Vi",          # Q -> Flash
        "Jarvan",      # E+Q -> Flash
        "Galio",       # E -> Flash
        "Ornn",        # E -> Flash
        "Rakan",       # W -> Flash
        "Pyke",        # E -> Flash (repositions Pyke, changing the stun path)
        "Rell",        # W (Crash Down) -> Flash
        "Riven",       # Q3 -> Flash
        "Sejuani",     # Q -> Flash
        "Shen",        # E -> Flash
        "Tryndamere",  # E -> Flash
        "Urgot",       # E -> Flash
        "Yone",        # Q3 -> Flash
        "Zac",         # E -> Flash
        "Camille"      # E2 (Wall Dive) -> Flash
    ]

    # The problem asks for the output separated by ",". The space after the comma improves readability.
    result = ", ".join(champions)
    print(result)

find_dash_flash_champions()