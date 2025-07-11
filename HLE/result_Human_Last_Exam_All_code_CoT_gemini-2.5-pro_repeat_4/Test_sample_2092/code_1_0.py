def solve():
    """
    This function identifies and prints the names of League of Legends champions
    from a predefined list who can perform the R-Flash combo.
    """
    # List of champions confirmed to be able to perform the R-Flash combo
    # based on game mechanics up to Season 14 (2024).
    # The combo involves casting R and then Flashing during the cast time
    # to reposition the ability's effect.
    r_flash_champions = [
        "Amumu",       # R cast time allows Flash to reposition the AoE root.
        "Azir",        # Can R and Flash behind enemies to shuffle them.
        "Braum",       # Can R and Flash to change the starting point of the fissure.
        "Cassiopeia",  # Can R and Flash to reposition the petrifying gaze cone.
        "Diana",       # Can R and Flash during the cast time to move the pull center.
        "Gnar",        # (Mega Gnar) Can R and Flash to change the knockback direction.
        "Hecarim",     # Can Flash during the R charge to change the fear location.
        "Lee Sin",     # Can R a target and Flash to change the kick direction.
        "Nami",        # Can R and Flash to have the tidal wave originate from the new spot.
        "Neeko",       # Can R and Flash just before landing to move the stun AoE.
        "Nilah",       # Can R and Flash during the cast time to move the pull center.
        "Orianna",     # With the ball on herself, can R and Flash to move the shockwave.
        "Qiyana",      # Can R and Flash during the cast time to reposition the knockback.
        "Rell",        # Can R and Flash to pull enemies from a surprise location.
        "Riven",       # Can use 2nd R (Wind Slash) and Flash to change the projectile origin.
        "Sejuani",     # Can R and Flash during the cast animation to throw from a new spot.
        "Seraphine",   # Can R and Flash to have the charm wave originate from a new spot.
        "Sona",        # Can R and Flash to reposition the stun AoE.
        "Thresh",      # Can R and Flash to reposition The Box before walls appear.
        "Xin Zhao",    # Can R and Flash to apply the knockback from the new location.
        "Yone"         # Can R and Flash during the channel to alter the dash's starting point/angle.
    ]

    # Print the final list, joined by commas.
    print(", ".join(r_flash_champions))

solve()
<<<Amumu, Azir, Braum, Cassiopeia, Diana, Gnar, Hecarim, Lee Sin, Nami, Neeko, Nilah, Orianna, Qiyana, Rell, Riven, Sejuani, Seraphine, Sona, Thresh, Xin Zhao, Yone>>>