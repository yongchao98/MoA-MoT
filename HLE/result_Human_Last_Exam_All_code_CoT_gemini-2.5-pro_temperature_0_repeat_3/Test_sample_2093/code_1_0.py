def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a buffered Dash-Flash combo.
    """
    # List of champions who can perform the Dash-Flash combo as per the specified criteria.
    # The combo involves using Flash during the dash animation to alter the ability's landing spot and effect.
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
        "Zac"
    ]

    # Join the list of names into a single string, separated by ", "
    result = ", ".join(champions)
    
    print(result)

find_dash_flash_champions()