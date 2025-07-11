def find_dash_flash_champions():
    """
    This function identifies champions from a predefined list who can perform a "Dash-Flash" combo.
    The combo is defined as using Flash during a dash ability's animation to buffer it,
    extending its range or changing its trajectory to land the ability's primary effect (damage/CC).
    """
    
    # List of champions who can perform the Dash-Flash combo as described.
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Ornn",
        "Rakan",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Camille"
    ]
    
    # Join the list of names into a single string, separated by ", "
    result = ", ".join(champions)
    
    # Print the final result
    print(result)

find_dash_flash_champions()