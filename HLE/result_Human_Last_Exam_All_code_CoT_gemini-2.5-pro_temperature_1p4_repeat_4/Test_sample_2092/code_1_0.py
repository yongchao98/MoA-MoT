def solve():
    """
    This function identifies and prints the names of League of Legends champions from a given list
    who can perform an R-Flash combo up until Season 14 (2024).
    """
    champions_who_can_r_flash = [
        "Amumu", "Aphelios", "Braum", "Cassiopeia", "Diana", "Fizz", "Gnar",
        "Graves", "Illaoi", "Irelia", "Janna", "Lee Sin", "Nami", "Neeko",
        "Nilah", "Qiyana", "Rell", "Renata Glasc", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Xin Zhao"
    ]
    
    result = ", ".join(champions_who_can_r_flash)
    print(result)

solve()