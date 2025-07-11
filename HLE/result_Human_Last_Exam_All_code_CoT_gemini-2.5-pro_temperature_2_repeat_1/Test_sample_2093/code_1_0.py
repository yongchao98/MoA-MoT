def solve():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform a 'Dash-Flash' combo. The names are printed as a single,
    comma-separated string.
    """
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
        "Camille"
    ]
    
    # The problem asks to print the answer, so we will use the print function.
    # The champion names should be separated by ",".
    result = ", ".join(champions)
    print(result)

solve()