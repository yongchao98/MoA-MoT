def solve():
    """
    This function identifies and prints the names of League of Legends champions
    who can perform a "Dash-Flash" combo with their specified abilities.
    The Dash-Flash combo involves using Flash during the dash animation to alter
    its endpoint and effect.
    """
    
    # List of champions that can perform the specified Dash-Flash combo.
    # Analysis of each champion:
    # 1. Gragas (E): Yes, can E and Flash during the dash to change the impact location.
    # 2. Vi (Q): Yes, can Flash during the Q dash to redirect the knock-up.
    # 3. Jarvan (E+Q): Yes, can Flash during the Q dash towards the flag to change the knock-up location.
    # 4. Fizz (E): No, Flash cannot be used during the damage-dealing part of the dash to redirect it.
    # 5. LeBlanc (W): No, W is an instant blink; no dash animation to Flash during.
    # 6. Galio (E): Yes, can Flash during the forward dash part to redirect the knock-up.
    # 7. Kayn (Q): No, the dash is too quick to buffer a Flash.
    # 8. Ornn (E): Yes, can E and Flash into terrain to create a surprise knock-up.
    # 9. Poppy (E): No, E is a targeted ability and cannot be redirected with Flash mid-dash.
    # 10. Rakan (W): Yes, can Flash just before landing to change the knock-up area.
    # 11. Pyke (E): No, Pyke's dash completes, and then Flash is used to redirect the returning phantom, not the dash itself.
    # 12. Rell (W): Yes, can Flash just before landing from her leap to change the knock-up area.
    # 13. Riven (Q3): Yes, can Flash during the third Q's hop to redirect the knock-up.
    # 14. Sejuani (Q): Yes, can Flash during the dash to change its trajectory.
    # 15. Shen (E): Yes, can E-Flash to extend the taunt range.
    # 16. Tristana (W): No, the jump's trajectory is locked once cast.
    # 17. Tryndamere (E): Yes, can E-Flash to extend the spin's damage range.
    # 18. Urgot (E): Yes, can Flash during his dash to change his final position and flip an enemy.
    # 19. Yone (Q3): Yes, can Flash during the third Q dash to redirect the knock-up.
    # 20. Zac (E): Yes, can Flash just before landing to change the knock-up location.
    # 21. Renekton (E): No, dash is too quick to buffer with Flash.
    # 22. Nidalee (W): No, dash is too quick to buffer.
    # 23. Fiora (Q): No, dash is too quick to buffer.
    # 24. Corki (W): No, the path is locked.
    # 25. Camille (E): No, the wall dive part of her E cannot be Flash-redirected.
    # 26. BelVeth (Q): No, dashes are too short and quick.

    champions_with_dash_flash = [
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
        "Zac"
    ]

    # Join the list of names with a comma and print the result.
    result = ",".join(champions_with_dash_flash)
    print(result)

solve()