def find_r_flash_champions():
    """
    Identifies and prints the names of League of Legends champions from a predefined list
    who can perform the R-Flash combo.
    """
    
    # List of champions provided by the user
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia",
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim",
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred",
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah",
        "Oriana", "Poppy", "Qiyana", "Rell", "Renata", "Riven (for the 2nd R)",
        "Sejuani", "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao",
        "Yone", "Zehri"
    ]
    
    # Champions who can perform the R-Flash combo as defined
    # R is cast, and Flash is used during the cast time to reposition the ability's effect.
    r_flash_champions = [
        "Amumu",        # R (Curse of the Sad Mummy) has a cast time, can Flash to move the AoE.
        "Azir",         # R (Emperor's Divide) can be cast, then Flash behind enemies to change shove direction.
        "Braum",        # R (Glacial Fissure) can be cast, then Flash to change the origin point of the fissure.
        "Cassiopeia",   # R (Petrifying Gaze) has a cast time, can Flash to change the cone's origin.
        "Gnar",         # Mega Gnar's R (GNAR!) can be cast, then Flash to reposition the knockback.
        "Graves",       # R (Collateral Damage) has a cast time, can Flash to fire from a new location.
        "Lee Sin",      # R (Dragon's Rage) can be cast, then Flash behind the target to change kick direction.
        "Nami",         # R (Tidal Wave) has a cast time, can Flash to change the wave's origin.
        "Neeko",        # R (Pop Blossom) can be channeled, then Flash right before it lands to move the AoE.
        "Nilah",        # R (Apotheosis) can be cast, then Flash during the spin to reposition the pull effect.
        "Qiyana",       # R (Supreme Display of Talent) has a cast time, can Flash to change the shockwave's angle.
        "Rell",         # R (Magnet Storm) can be cast, then Flash to pull enemies from a new position.
        "Sejuani",      # R (Glacial Prison) has a cast time, can Flash to throw the bola from a new spot.
        "Seraphine",    # R (Encore) has a cast time, can Flash to change the projectile's origin.
        "Skarner",      # Pre-rework R (Impale) allowed Flash during the suppression to drag the target.
        "Sona",         # R (Crescendo) has a cast time, can Flash to move the stun area.
        "Xin Zhao",     # R (Crescent Guard) can be cast, then Flash during the sweep animation to change the knockback center.
        "Yone"          # R (Fate Sealed) has a channel, can Flash during it to reposition the entire dash.
    ]

    print(", ".join(r_flash_champions))

find_r_flash_champions()