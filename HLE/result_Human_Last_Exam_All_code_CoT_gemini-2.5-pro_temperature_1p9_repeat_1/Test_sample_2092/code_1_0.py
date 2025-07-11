def solve_r_flash_champions():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform an R-Flash combo as of Season 14 (2024).
    The R-Flash combo is defined as casting R then using Flash during the cast time
    to alter the ability's effect for optimal impact.
    """
    
    # List of champions that can perform the R-Flash combo based on the specific definition.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Gnar",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Rell",
        "Riven (for the 2nd R)",
        "Sejuani",
        "Seraphine",
        "Sona",
        "Xin Zhao",
        "Yone",
        "Zeri" # Correcting the typo from 'Zehri'
    ]
    
    # Join the names with a comma and space for clean output
    result = ", ".join(r_flash_champions)
    
    print(result)

solve_r_flash_champions()