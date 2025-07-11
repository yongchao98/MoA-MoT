def solve():
    """
    This function identifies and prints the names of League of Legends champions
    who can perform a specific buffered "Dash-Flash" combo.
    """
    champions = [
        "Gragas", "Vi", "Jarvan", "Galio", "Kayn",
        "Ornn", "Rakan", "Pyke", "Rell", "Riven",
        "Sejuani", "Shen", "Tryndamere", "Urgot",
        "Yone", "Fiora", "Camille"
    ]
    
    # The user requested the names to be separated by a comma.
    print(",".join(champions))

solve()