def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a Dash-Flash combo.
    """
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Orn",
        "Rakan",
        "Pyke",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tristana",
        "Tryndamere",
        "Urgot",
        "Zac",
        "Renekton"
    ]
    
    # Joining the list of champion names into a single string, separated by commas.
    result = ", ".join(champions)
    
    print(result)

find_dash_flash_champions()