def find_dash_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform a Dash-Flash combo.
    """
    champions = {
        "Gragas": True,
        "Vi": True,
        "Jarvan": True,
        "Fizz": False,
        "LeBlanc": False,
        "Galio": True,
        "Kayn": False,
        "Ornn": True,
        "Poppy": False,
        "Rakan": True,
        "Pyke": False,
        "Rell": True,
        "Riven": True,
        "Sejuani": True,
        "Shen": True,
        "Tristana": False,
        "Tryndamere": True,
        "Urgot": True,
        "Yone": True,
        "Zac": True,
        "Renekton": False,
        "Nidalee": False,
        "Fiora": False,
        "Corki": False,
        "Camille": True,
        "BelVeth": False,
    }

    valid_champions = [name for name, can_combo in champions.items() if can_combo]
    
    print(",".join(valid_champions))

find_dash_flash_champions()