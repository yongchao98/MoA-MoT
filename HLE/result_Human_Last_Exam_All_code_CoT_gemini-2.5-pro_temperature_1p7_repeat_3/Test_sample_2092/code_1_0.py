def find_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    The analysis is based on game mechanics knowledge up to Season 14 (2024).
    """
    champions = [
        "Amumu", "Aphelios", "Azir", "Blitzcrank", "Braum", "Cassiopeia",
        "Diana", "Evelynn", "Fizz", "Gnar", "Graves", "Gwen", "Hecarim",
        "Illaoi", "Irelia", "Janna", "Jarvan", "Jax", "Kassadin", "Kindred",
        "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus", "Neeko", "Nilah",
        "Orianna", "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone",
        "Zehri" # Note: This seems to be a typo for Zeri
    ]

    # Based on game knowledge, these champions can perform the R-Flash combo
    # where Flash is used during the R cast/channel time to change the
    # ability's effective origin or area.
    r_flash_champions = [
        "Amumu",        # R has a cast time, Flash moves the center of the explosion.
        "Braum",        # R has a cast time, Flash moves the origin of the fissure.
        "Cassiopeia",   # R has a cast time, Flash moves the origin of the cone.
        "Diana",        # R has a short delay, Flash can reposition the center of the pull.
        "Gnar",         # Mega Gnar's R has a cast time, Flash repositions him for the shove.
        "Graves",       # R has a cast time, Flash can cancel the self-knockback.
        "Janna",        # R is a channel, Flash can reposition the knockback effect.
        "Lee Sin",      # R has a cast time, Flash repositions Lee for the kick direction (Insec).
        "Nami",         # R has a cast time, Flash moves the origin of the wave.
        "Neeko",        # R is a channel, Flash moves the center of the explosion right before it lands.
        "Nilah",        # R has a cast time, Flash can reposition the center of the pull.
        "Qiyana",       # R has a cast time, Flash can change the angle of the shockwave.
        "Rell",         # R is a channel, Flash can pull enemies from a new position.
        "Renata",       # R has a cast time, Flash can move the origin of the projectile wave.
        "Riven",        # The second R (Wind Slash) has a cast time, Flash moves the projectile origin.
        "Sejuani",      # R has a cast time, Flash moves the origin of the thrown bola.
        "Seraphine",    # R has a cast time, Flash moves the origin of the projectile.
        "Skarner",      # Pre-rework R allowed Flash during suppression to drag enemies further.
        "Sona",         # R has a cast time, Flash moves the origin of the stun area.
        "Xin Zhao",     # R has a cast time, Flash can reposition the center of the sweep/knockback.
        "Yone"          # R has a channel, Flash can change the start position of the dash.
    ]

    print(", ".join(r_flash_champions))

find_r_flash_champions()
<<<Amumu, Braum, Cassiopeia, Diana, Gnar, Graves, Janna, Lee Sin, Nami, Neeko, Nilah, Qiyana, Rell, Renata, Riven, Sejuani, Seraphine, Skarner, Sona, Xin Zhao, Yone>>>