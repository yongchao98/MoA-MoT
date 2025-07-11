def find_dash_flash_champions():
    """
    Identifies and prints the names of champions from a predefined list who can
    perform an effective Dash-Flash combo to land their ability's effect.
    """
    # List of champions who can perform the specified Dash-Flash combo
    champions = [
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
        "Renekton",
        "Fiora",
        "Camille"
    ]

    # Join the list into a single string, separated by ", "
    result = ", ".join(champions)

    # Print the final result
    print(result)

find_dash_flash_champions()
