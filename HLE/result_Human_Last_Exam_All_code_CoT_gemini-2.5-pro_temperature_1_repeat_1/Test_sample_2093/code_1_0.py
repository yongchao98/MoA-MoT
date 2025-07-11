def find_dash_flash_champions():
    """
    This function identifies champions from a predefined list who can perform a "Dash-Flash" combo.
    The combo is defined as using Flash during a dash ability's animation to redirect or extend its effect (damage/CC),
    not just for repositioning.
    """
    
    # List of all champions and their abilities to check
    # 1. Gragas (E) - Yes
    # 2. Vi (Q) - Yes
    # 3. Jarvan (E+Q) - Yes
    # 4. Fizz (E) - No
    # 5. LeBlanc (W) - No
    # 6. Galio (E) - Yes
    # 7. Kayn (Q) - No
    # 8. Ornn (E) - Yes
    # 9. Poppy (E) - No (targeted)
    # 10. Rakan (W) - Yes
    # 11. Pyke (E) - Yes
    # 12. Rell (W) - Yes
    # 13. Riven (Q3) - Yes
    # 14. Sejuani (Q) - Yes
    # 15. Shen (E) - Yes
    # 16. Tristana (W) - No
    # 17. Tryndamere (E) - Yes
    # 18. Urgot (E) - Yes
    # 19. Yone (Q3) - Yes
    # 20. Zac (E) - Yes
    # 21. Renekton (E) - No
    # 22. Nidalee (W) - No
    # 23. Fiora (Q) - No
    # 24. Corki (W) - No
    # 25. Camille (E) - Yes
    # 26. BelVeth (Q) - No

    champions_who_can_dash_flash = [
        "Gragas",
        "Vi",
        "Jarvan",
        "Galio",
        "Ornn",
        "Rakan",
        "Pyke",
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

    # Join the list of names with ", "
    result = ", ".join(champions_who_can_dash_flash)
    
    # Print the final result
    print(result)

find_dash_flash_champions()