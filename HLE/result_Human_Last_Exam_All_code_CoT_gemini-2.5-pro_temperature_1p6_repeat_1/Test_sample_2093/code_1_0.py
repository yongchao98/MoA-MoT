def find_dash_flash_champions():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform a 'Dash-Flash' combo to buffer their ability's effect.
    """
    
    # List of champions and abilities provided by the user.
    # The analysis is done based on game knowledge.
    # 1. Gragas (E): Yes, E-Flash changes the hitbox location.
    # 2. Vi (Q): Yes, Q-Flash changes the impact location.
    # 3. Jarvan (E+Q): Yes, E-Q-Flash extends the knock-up range.
    # 4. Fizz (E): No, cannot flash during E to change its impact zone.
    # 5. LeBlanc (W): No, W is a blink, not a dash; no travel time to Flash during.
    # 6. Galio (E): Yes, E-Flash changes the knock-up location.
    # 7. Kayn (Q): No, dash is too fast and cannot be buffered with Flash.
    # 8. Ornn (E): Yes, E-Flash can redirect him into a wall for a knock-up.
    # 9. Poppy (E): No, it's a targeted dash, cannot be redirected with Flash.
    # 10. Rakan (W): Yes, W-Flash changes the knock-up location.
    # 11. Pyke (E): Yes, E-Flash extends the stun range of the phantom.
    # 12. Rell (W): Yes, W-Flash changes the landing/knock-up spot.
    # 13. Riven (Q3): Yes, Q3-Flash changes the landing/knock-up spot.
    # 14. Sejuani (Q): Yes, Q-Flash changes the dash path and impact point.
    # 15. Shen (E): Yes, E-Flash changes the dash path and taunt area.
    # 16. Tristana (W): Yes, W-Flash changes the landing spot and damage/slow area.
    # 17. Tryndamere (E): Yes, E-Flash extends the spin's damage range.
    # 18. Urgot (E): Yes, E-Flash changes the collision point for the stun/flip.
    # 19. Yone (Q3): Yes, Q3-Flash changes the dash path and knock-up area.
    # 20. Zac (E): Yes, can Flash mid-flight to change landing and knock-up spot.
    # 21. Renekton (E): No, dash is too fast and cannot be buffered.
    # 22. Nidalee (W): No, cougar pounce cannot be buffered with Flash.
    # 23. Fiora (Q): No, cannot Flash during the lunge to redirect its effect.
    # 24. Corki (W): Yes, can Flash during Valkyrie to change path and damage trail.
    # 25. Camille (E): Yes, E2-Flash extends the stun range.
    # 26. BelVeth (Q): No, dashes are too fast and cannot be buffered.

    champions_who_can_combo = [
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
        "Tristana",
        "Tryndamere",
        "Urgot",
        "Yone",
        "Zac",
        "Corki",
        "Camille"
    ]
    
    result_string = ", ".join(champions_who_can_combo)
    print(result_string)

find_dash_flash_champions()