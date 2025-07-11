def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a true "Dash-Flash" combo.
    The combo involves using Flash during a dash ability to change the
    location of the ability's effect (e.g., damage or crowd control).
    """
    
    # List of champions who can perform the specified Dash-Flash combo.
    # Analysis:
    # Gragas (E), Vi (Q), Jarvan (E+Q), Galio (E), Ornn (E), Rakan (W), Rell (W),
    # Riven (Q3), Sejuani (Q), Shen (E), Tryndamere (E), Urgot (E), Yone (Q3),
    # Zac (E), and Camille (E2) can all buffer Flash during their dash
    # to change the ability's impact zone.
    # Others like LeBlanc, Tristana, or Corki can Flash during their dash,
    # but the ability's damage lands at the original location, only repositioning the champion.
    # Fizz, Poppy, Pyke, and Fiora's abilities do not work with this specific mechanic.
    
    champions = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Ornn",
        "Rakan",
        "Rell",
        "Riven",
        "Sejuani",
        "Shen",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Camille"
    ]
    
    # Join the list into a single string, separated by commas.
    result = ",".join(champions)
    
    print(result)

find_dash_flash_champions()