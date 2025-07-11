def find_dash_flash_champions():
    """
    This function identifies and prints the names of champions who can perform a Dash-Flash combo.
    The list is predetermined based on in-game mechanics.
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
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Camille"
    ]
    
    result = ",".join(champions)
    print(result)

find_dash_flash_champions()