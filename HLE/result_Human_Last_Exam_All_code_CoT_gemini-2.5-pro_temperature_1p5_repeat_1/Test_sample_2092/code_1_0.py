def find_r_flash_champions():
    """
    This function identifies League of Legends champions from a predefined list
    who can perform the R-Flash combo up until Season 14 (2024).

    The R-Flash combo is defined as casting the ultimate (R) and then using Flash
    during the cast time to buffer the ability, causing it to resolve from the
    new location. This is for a more optimal landing of damage or crowd-control.
    """

    # The list of champions provided by the user.
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana",
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi",
        "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin",
        "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana",
        "Poppy", "Qiyana", "Rell", "Renata", "Riven (for the 2nd R)", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # Champions who can perform the R-Flash combo based on ability mechanics.
    # Note: Riven is simplified to "Riven" for the final list.
    r_flash_capable = [
        "Amumu",      # R has a cast time, Flash repositions the AoE.
        "Azir",       # Classic 'Shurima Shuffle', R repositions the soldier wall.
        "Braum",      # R has a cast time, Flash repositions the fissure's origin.
        "Cassiopeia", # R has a cast time, Flash repositions the cone.
        "Gnar",       # Mega Gnar's R can be Flashed to change the shove angle/position.
        "Lee Sin",    # The iconic 'Insec', R a target and Flash behind them to kick them differently.
        "Neeko",      # Can R and Flash before landing to reposition the AoE stun.
        "Poppy",      # Tap R has a short animation that can be buffered with Flash.
        "Qiyana",     # R has a cast time, Flash repositions Qiyana, altering the shockwave.
        "Riven",      # The second activation of R (Wind Slash) can be R-Flashed to redirect.
        "Sejuani",    # R has a cast time, Flash repositions the projectile's origin point.
        "Seraphine",  # R has a cast time, Flash repositions the charm projectile's origin.
        "Skarner",    # Both pre and post-rework R can be buffered with Flash.
        "Sona",       # R has a cast time, Flash repositions the AoE stun.
        "Xin Zhao",   # R has a cast time for its sweep, Flash repositions the knockback.
        "Yone"        # R has a wind-up, Flash repositions the start of his dash.
    ]

    # The prompt included "Riven (for the 2nd R)". We'll use just "Riven" in the output.
    final_list = []
    for champ in r_flash_capable:
        if "Riven" in champ:
            final_list.append("Riven")
        else:
            final_list.append(champ)

    # Print the final list, with each champion name separated by a comma.
    print(", ".join(final_list))

find_r_flash_champions()