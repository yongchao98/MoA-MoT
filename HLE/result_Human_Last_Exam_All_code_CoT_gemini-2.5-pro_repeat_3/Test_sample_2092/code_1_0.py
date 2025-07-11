def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    """
    # The full list of champions provided by the user.
    # Note: Some names have been corrected for accuracy (e.g., Jarvan IV, Zeri).
    champion_pool = [
        "Amumu", "Aphelios", "Azir", "Blitzcrank", "Braum", "Cassiopeia",
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim",
        "Illaoi", "Irelia", "Janna", "Jarvan IV", "Jax", "Kassadin",
        "Kindred", "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus",
        "Neeko", "Nilah", "Orianna", "Poppy", "Qiyana", "Rell", "Renata Glasc",
        "Riven", "Sejuani", "Seraphine", "Skarner", "Sona", "Thresh", "Vi",
        "Xin Zhao", "Yone", "Zeri"
    ]

    # List of champions who can perform the R-Flash combo as described.
    # The condition is that R is cast first, then Flash is used during the
    # cast time/buffer to alter the ability's point of origin.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Diana",
        "Evelynn",
        "Gnar",
        "Hecarim",
        "Illaoi",
        "Jarvan IV",
        "Lee Sin",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Xin Zhao",
        "Yone"
    ]

    # Print the final list, with names separated by a comma and a space.
    print(",".join(r_flash_champions))

solve()
<<<Amumu,Azir,Braum,Cassiopeia,Diana,Evelynn,Gnar,Hecarim,Illaoi,Jarvan IV,Lee Sin,Neeko,Nilah,Qiyana,Rell,Riven,Sejuani,Seraphine,Skarner,Sona,Xin Zhao,Yone>>>