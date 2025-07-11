def find_r_flash_champions():
    """
    This function identifies champions from a predefined list who can perform an R-Flash combo.

    An R-Flash is a technique where a player casts their ultimate ability (R) and
    then immediately uses Flash during the cast time to reposition the ability's effect.
    """

    # The full list of champions provided by the user.
    all_champions = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana",
        "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi", "Irelia",
        "Janna", "Jarvan", "Jax", "Kassadin", "Kindred", "Lee Sin", "Malphite",
        "Milio", "Nami", "Nautilus", "Neeko", "Nilah", "Oriana", "Poppy",
        "Qiyana", "Rell", "Renata", "Riven", "Sejuani", "Seraphine", "Skarner",
        "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # Champions who can perform the R-Flash combo based on game mechanics up to Season 14.
    # The logic for inclusion is having an ultimate with a cast time where Flash can
    # be used to change the ability's origin point or effective area.
    r_flash_capable = [
        "Amumu",       # R has a cast time; Flash moves the AOE's center.
        "Azir",        # R can be Flashed during cast to push enemies from a new angle.
        "Braum",       # R has a cast time; Flash moves the fissure's starting point.
        "Cassiopeia",  # R has a cast time; Flash repositions the cone of effect.
        "Gnar",        # Mega Gnar R can be Flashed to reposition the shove.
        "Graves",      # R has a cast time; Flash changes the projectile's origin.
        "Hecarim",     # R fear originates from his landing position, which can be changed with Flash.
        "Jarvan",      # R leap can be redirected with Flash before landing.
        "Lee Sin",     # R can be Flashed during the kick animation to change the kick's direction.
        "Neeko",       # R channel allows for Flash to reposition the final crash.
        "Orianna",     # R shockwave can be moved if the ball is on her and she Flashes.
        "Qiyana",      # R has a cast time; Flash repositions the initial shockwave.
        "Rell",        # R pull can be repositioned with Flash.
        "Riven",       # Second R (Wind Slash) has a cast time; Flash changes projectile origin.
        "Sejuani",     # R has a cast time; Flash changes the projectile's origin.
        "Seraphine",   # R has a cast time; Flash changes the projectile's origin.
        "Sona",        # R has a cast time; Flash repositions the stun area.
        "Vi",          # R can be Flashed just before the dash to alter the path and knockup angle.
        "Xin Zhao",    # R has a cast time; Flash repositions the knockback sweep.
        "Yone"         # R has a cast time; Flash repositions the start of his dash.
    ]

    # Print the final list as a comma-separated string.
    print(",".join(r_flash_capable))

find_r_flash_champions()
<<<Amumu,Azir,Braum,Cassiopeia,Gnar,Graves,Hecarim,Jarvan,Lee Sin,Neeko,Orianna,Qiyana,Rell,Riven,Sejuani,Seraphine,Sona,Vi,Xin Zhao,Yone>>>