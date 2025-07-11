def solve():
    """
    This function identifies and prints the names of League of Legends champions
    who can perform a 'Dash-Flash' combo as defined in the prompt.
    """
    champions_who_can_dash_flash = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Orn",
        "Rakan",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tristana",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Renekton",
        "Corki",
        "Camille"
    ]
    
    # Joining the list of names with a comma separator for the final output
    result = ",".join(champions_who_can_dash_flash)
    
    print(result)

solve()