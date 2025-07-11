def solve_r_flash_champions():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    The R-Flash combo is defined as casting the ultimate ability (R) and then
    using Flash during the cast time to redirect or reposition the ability's effect.
    """

    # The list is based on analyzing each champion's ultimate ability (R) up to Season 14 (2024).
    # Champions are included if their R has a cast time/channel during which Flash can
    # be used to change the ability's origin or area of effect.
    r_flash_champions = [
        "Amumu",      # R has a cast time; Flash moves the center of the curse.
        "Azir",       # R has a cast time; Flash changes the starting position of the soldier wall.
        "Braum",      # R has a cast time; Flash changes the origin point of the fissure.
        "Cassiopeia", # R has a cast time; Flash moves the cone of her gaze.
        "Diana",      # R has a charge time; Flash moves the center of the pull/damage.
        "Gnar",       # Mega Gnar's R has a cast time; Flash repositions the shove.
        "Lee Sin",    # R has a cast animation; Flash repositions Lee to kick the target in a new direction.
        "Neeko",      # R has a channel; Flash repositions the center of the final stun.
        "Nilah",      # R has a cast time; Flash repositions the center of the pull.
        "Poppy",      # Charged R has a channel; Flash changes the angle of the swing.
        "Qiyana",     # R has a cast time; Flash changes the origin point of the shockwave.
        "Rell",       # R has an initial cast; Flash can be used to drag pulled enemies.
        "Riven",      # The second activation of R (Wind Slash) has a cast time; Flash changes projectile origins.
        "Sejuani",    # R has a throw animation; Flash changes the origin point of the projectile.
        "Seraphine",  # R has a cast time; Flash changes the origin point of the charm projectile.
        "Skarner",    # R has a windup; Flash can change the skillshot's origin (applies to pre and post-rework).
        "Sona",       # R has a very short (0.25s) cast time; a frame-perfect Flash can move the Crescendo's origin.
        "Xin Zhao",   # R has a cast time; Flash can change the center of the knockback sweep.
        "Yone"        # R has a channel before the dash; Flash can change the starting point of the charge.
    ]

    print(",".join(r_flash_champions))

solve_r_flash_champions()