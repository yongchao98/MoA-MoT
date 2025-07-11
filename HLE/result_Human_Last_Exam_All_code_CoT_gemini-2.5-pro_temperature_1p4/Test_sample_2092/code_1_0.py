def find_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    The R-Flash technique involves casting the ultimate ability (R) and then
    using Flash during the cast time to change the ability's outcome.
    """
    
    # List of champions who can buffer their Ultimate with Flash
    # based on game mechanics up to Season 14 (2024).
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Hecarim",
        "Lee Sin",
        "Neeko",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Thresh",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]
    
    # Print the champion names separated by a comma
    print(",".join(r_flash_champions))

find_r_flash_champions()
<<<Amumu,Azir,Braum,Cassiopeia,Gnar,Hecarim,Lee Sin,Neeko,Rell,Riven,Sejuani,Seraphine,Skarner,Sona,Thresh,Vi,Xin Zhao,Yone>>>