def find_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform an R-Flash combo.
    """
    champion_list = [
        "Amumu", "Aphelios", "Azir", "Blitzcrank", "Braum", "Cassiopeia",
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim",
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred",
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah",
        "Orianna", "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # List of champions who can perform a true R-Flash combo
    # where Flash is used during the R cast/channel to alter its origin/effect.
    r_flash_champions = [
        "Amumu",
        "Azir",
        "Braum",
        "Cassiopeia",
        "Evelynn",
        "Gnar",
        "Graves",
        "Hecarim",
        "Jarvan",
        "Lee Sin",
        "Neeko",
        "Orianna",
        "Qiyana",
        "Rell",
        "Riven",
        "Sejuani",
        "Seraphine",
        "Skarner",
        "Sona",
        "Vi",
        "Xin Zhao",
        "Yone"
    ]

    print(",".join(r_flash_champions))

find_r_flash_champions()
<<<Amumu,Azir,Braum,Cassiopeia,Evelynn,Gnar,Graves,Hecarim,Jarvan,Lee Sin,Neeko,Orianna,Qiyana,Rell,Riven,Sejuani,Seraphine,Skarner,Sona,Vi,Xin Zhao,Yone>>>