def find_dash_flash_champions():
    """
    This function identifies and prints the names of champions from a predefined list
    who can perform a 'Dash-Flash' combo.
    """
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Orn",
        "Pyke",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Camille"
    ]
    
    result = ",".join(champions)
    print(result)

find_dash_flash_champions()