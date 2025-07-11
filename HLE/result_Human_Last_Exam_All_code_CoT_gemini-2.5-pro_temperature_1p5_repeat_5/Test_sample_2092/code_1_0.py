def solve():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform the R-Flash combo up until season 14 (2024).
    The R-Flash combo is defined as using Flash during the ultimate's cast time
    to optimally change its area of effect, origin, or direction.
    """
    
    # List of all champions provided in the problem
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana", 
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi", "Irelia", 
        "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin", "Malphite", 
        "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana", "Poppy", "Qiyana", 
        "Rell", "Renata", "Riven (for the 2nd R)", "Sejuani", "Seraphine", "Skarner", 
        "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # Champions who can perform the R-Flash combo as described.
    # Note: 'Zehri' from the original list is assumed to be Zeri.
    r_flash_champions = [
        "Amumu",
        "Braum",
        "Cassiopeia",
        "Gnar",
        "Graves",
        "Lee Sin",
        "Nami",
        "Neeko",
        "Nilah",
        "Qiyana",
        "Rell",
        "Sejuani",
        "Seraphine",
        "Sona",
        "Xin Zhao",
        "Yone",
        "Zehri"
    ]

    # Print the final list of champion names, separated by a comma.
    print(", ".join(r_flash_champions))

solve()