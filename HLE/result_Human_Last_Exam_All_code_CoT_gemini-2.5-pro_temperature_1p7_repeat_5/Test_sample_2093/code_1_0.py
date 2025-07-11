def find_dash_flash_champions():
    """
    This function identifies champions from a predefined list who can perform a "Dash-Flash" combo.

    The "Dash-Flash" combo is defined as using Flash during a dash ability's animation
    to change the location where the ability's primary effect (damage/CC) is applied.
    """

    # List of champions and their abilities to be checked.
    # Note: This analysis is based on established game mechanics.
    # 1. Gragas (E): Yes, classic example.
    # 2. Vi (Q): Yes, can flash during the Q dash.
    # 3. Jarvan (E+Q): Yes, can flash during the Q travel to the flag.
    # 4. Fizz (E): No, you can't flash during the dash portion to change its impact zone.
    # 5. LeBlanc (W): No, W is a blink, not a dash with travel time to buffer Flash in.
    # 6. Galio (E): Yes, can flash during the forward punch.
    # 7. Kayn (Q): No, the dash is too quick and you cannot buffer Flash to change the spin location.
    # 8. Ornn (E): Yes, can flash during the charge.
    # 9. Poppy (E): No, E is a targeted dash, cannot be redirected with Flash.
    # 10. Rakan (W): Yes, can flash during the leap to change the knock-up location.
    # 11. Pyke (E): No, the dash itself doesn't do CC/damage. Flash repositions Pyke for the phantom's return path, a different mechanic.
    # 12. Rell (W - Crash Down): Yes, can flash mid-air to change the landing spot.
    # 13. Riven (Q3): Yes, can flash during the third Q's hop.
    # 14. Sejuani (Q): Yes, can flash during the dash.
    # 15. Shen (E): Yes, can flash during the dash to taunt at the new location.
    # 16. Tristana (W): No, jump destination is locked on cast.
    # 17. Tryndamere (E): Yes, can flash mid-spin to apply damage elsewhere.
    # 18. Urgot (E): Yes, can flash during the dash.
    # 19. Yone (Q3): Yes, can flash during the dash to reposition the knock-up.
    # 20. Zac (E): Yes, can flash mid-flight.
    # 21. Renekton (E): No, you can't buffer flash during the dash itself.
    # 22. Nidalee (W): No, the pounce cannot be redirected mid-air with Flash.
    # 23. Fiora (Q): No, lunge path is locked on cast.
    # 24. Corki (W): No, Valkyrie path is locked on cast.
    # 25. Camille (E2): Yes, can flash during the second E (Wall Dive).
    # 26. BelVeth (Q): No, dashes are too fast to buffer a Flash combo.

    champions_who_can_dash_flash = [
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

    print(",".join(champions_who_can_dash_flash))

find_dash_flash_champions()