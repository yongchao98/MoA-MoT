def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform an R-Flash combo.
    The analysis is based on game mechanics knowledge up to season 14.
    """
    champions = [
        # Champions who CAN R-Flash
        "Amumu", "Braum", "Cassiopeia", "Gnar", "Graves", "Jarvan", "Lee Sin", 
        "Neeko", "Nilah", "Orianna", "Poppy", "Qiyana", "Riven", "Sejuani", 
        "Seraphine", "Skarner", "Sona", "Thresh", "Xin Zhao", "Yone"
    ]
    
    # Print the final list, with names separated by a comma.
    print(",".join(champions))

solve_r_flash_champions()