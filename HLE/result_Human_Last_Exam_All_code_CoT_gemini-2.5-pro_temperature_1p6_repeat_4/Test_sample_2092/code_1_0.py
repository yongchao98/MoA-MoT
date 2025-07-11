def find_r_flash_champions():
    """
    Identifies and prints the names of champions from a predefined list who can perform an R-Flash combo.
    The list is based on in-game mechanics as of early 2024.
    """
    champions = [
        "Amumu", "Azir", "Braum", "Cassiopeia", "Diana", "Evelynn", "Gnar",
        "Hecarim", "Janna", "Jarvan", "Lee Sin", "Neeko", "Nilah", "Orianna",
        "Poppy", "Qiyana", "Rell", "Riven", "Sejuani", "Seraphine",
        "Skarner", "Sona", "Xin Zhao", "Yone"
    ]
    
    result = ", ".join(champions)
    print(result)

find_r_flash_champions()