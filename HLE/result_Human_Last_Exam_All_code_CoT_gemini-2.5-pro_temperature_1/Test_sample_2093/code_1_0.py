def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a Dash-Flash combo.
    """
    
    # List of champions who can buffer Flash during their dash ability
    # to alter the ability's outcome.
    champions_who_can_dash_flash = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Ornn",
        "Rakan",
        "Pyke",
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
    
    # Join the list of names with ", " and print the result.
    result = ", ".join(champions_who_can_dash_flash)
    print(result)

find_dash_flash_champions()
<<<Gragas, Vi, Jarvan, Galio, Ornn, Rakan, Pyke, Rell, Riven, Sejuani, Shen, Tryndamere, Urgot, Yone, Zac, Camille>>>