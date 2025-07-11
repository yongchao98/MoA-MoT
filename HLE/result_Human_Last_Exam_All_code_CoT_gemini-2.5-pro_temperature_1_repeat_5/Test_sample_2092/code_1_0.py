def find_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions from a predefined list
    who can perform the R-Flash combo as of Season 14 (2024).
    """
    
    # The full list of champions provided in the prompt.
    champion_pool = [
        "Amumu", "Aphelios", "Azir", "Blitzrank", "Braum", "Cassiopeia", "Diana", "Evelynn",
        "Fizz", "Gnar", "Graves", "Gwen", "Hecarim", "Illaoi", "Irelia", "Janna", "Jarvan",
        "Jax", "Kassadin", "Kindred", "Lee Sin", "Malphite", "Milio", "Nami", "Nautilus",
        "Neeko", "Nilah", "Oriana", "Poppy", "Qiyana", "Rell", "Renata", "Riven", "Sejuani",
        "Seraphine", "Skarner", "Sona", "Thresh", "Vi", "Xin Zhao", "Yone", "Zehri"
    ]

    # Champions confirmed to be able to R-Flash.
    # The R ability has a cast time or channel that can be buffered by Flash,
    # causing the ability to originate from the new (Flashed) location.
    r_flash_champs = [
        "Amumu",      # R has a cast time, Flash relocates the AoE.
        "Braum",      # R has a cast time, Flash relocates the fissure's origin.
        "Cassiopeia", # R has a cast time, Flash relocates the cone's origin.
        "Diana",      # R has a channel, Flash relocates the pull/damage AoE.
        "Gnar",       # Mega Gnar R has a cast time, Flash can change the knockback position.
        "Graves",     # R has a cast time, Flash makes the projectile fire from the new location.
        "Lee Sin",    # R can be cast on a target, followed by a Flash to change the kick direction.
        "Nami",       # R has a cast time, Flash relocates the wave's origin.
        "Neeko",      # R has a channel, Flash right before it lands relocates the stun AoE.
        "Poppy",      # Tap R has a cast time, Flash relocates the knock-up.
        "Qiyana",     # R has a cast time, Flash relocates the shockwave's origin.
        "Rell",       # R is a channel, Flash can reposition the pulling effect.
        "Riven",      # R2 (Wind Slash) has a cast time, Flash makes projectiles fire from the new location.
        "Sejuani",    # R has a cast time, Flash makes the bola throw from the new location.
        "Seraphine",  # R has a cast time, Flash makes the projectile fire from the new location.
        "Skarner",    # The 2024 reworked R has a cast time, Flash can change the skillshot's origin.
        "Sona",       # R has a cast time, Flash relocates the stun AoE.
        "Xin Zhao",   # R has a cast time, Flash relocates the sweep AoE.
        "Yone"        # R has a channel, Flash during it relocates the entire dash and pull effect.
    ]

    print(",".join(r_flash_champs))

find_r_flash_champions()
<<<Amumu,Braum,Cassiopeia,Diana,Gnar,Graves,Lee Sin,Nami,Neeko,Poppy,Qiyana,Rell,Riven,Sejuani,Seraphine,Skarner,Sona,Xin Zhao,Yone>>>