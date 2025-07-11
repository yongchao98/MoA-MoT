def find_r_flash_champions():
    """
    This function identifies and prints the champions from a predefined list
    who can perform a buffered R-Flash combo in League of Legends.
    The R-Flash combo is defined as using Flash during the ultimate's cast
    time to alter its point of origin or effect area for optimal impact.
    """
    # List of champions (and their original numbers) who can perform the R-Flash combo
    # based on game mechanics up to Season 14 (2024).
    r_flash_champions = [
        (1, "Amumu"),
        (3, "Azir"),
        (5, "Braum"),
        (6, "Cassiopeia"),
        (7, "Diana"),
        (10, "Gnar"),
        (11, "Graves"),
        (21, "Lee Sin"),
        (23, "Milio"),
        (26, "Neeko"),
        (27, "Nilah"),
        (28, "Orianna"),
        (30, "Qiyana"),
        (31, "Rell"),
        (33, "Riven"),
        (34, "Sejuani"),
        (35, "Seraphine"),
        (36, "Skarner"),
        (37, "Sona"),
        (40, "Xin Zhao"),
        (41, "Yone"),
        (42, "Zehri")
    ]

    # Format the list to include the number and name for each champion
    formatted_list = [f"{num}. {name}" for num, name in r_flash_champions]

    # Join the formatted list into a single comma-separated string
    final_output_string = ", ".join(formatted_list)

    print(final_output_string)

find_r_flash_champions()