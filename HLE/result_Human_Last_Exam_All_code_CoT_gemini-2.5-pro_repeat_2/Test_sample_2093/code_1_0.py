def find_dash_flash_champions():
    """
    This function identifies champions from a predefined list who can perform
    a "Dash-Flash" combo, where Flash is used during the dash animation
    to alter the ability's effective landing spot for damage or crowd control.
    """
    champions = {
        1: "Gragas",
        2: "Vi",
        3: "Jarvan",
        4: "Fizz",
        5: "LeBlanc",
        6: "Galio",
        7: "Kayn",
        8: "Orn",
        9: "Poppy",
        10: "Rakan",
        11: "Pyke",
        12: "Rell",
        13: "Riven",
        14: "Sejuani",
        15: "Shen",
        16: "Tristana",
        17: "Tryndamere",
        18: "Urgot",
        19: "Yone",
        20: "Zac",
        21: "Renekton",
        22: "Nidalee",
        23: "Fiora",
        24: "Corki",
        25: "Camille",
        26: "BelVeth"
    }

    # List of champion numbers that can perform the Dash-Flash combo
    capable_champions_indices = [1, 2, 3, 6, 8, 10, 11, 12, 13, 14, 15, 17, 18, 20, 25]

    # Retrieve the names of the capable champions
    capable_champion_names = [champions[i] for i in capable_champions_indices]

    # Print the result
    print(",".join(capable_champion_names))

find_dash_flash_champions()