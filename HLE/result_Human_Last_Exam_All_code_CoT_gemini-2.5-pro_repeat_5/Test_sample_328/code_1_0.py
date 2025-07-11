def estimate_mating_age():
    """
    Provides the best average estimate of the male's age since eclosion
    for Tridactylophagus tartari and Periplaneta americana based on their life cycles.
    """
    
    # Species 1: Tridactylophagus tartari (a short-lived wasp)
    # Males emerge and mate soon after. With a lifespan of only a few days,
    # the average age of a mating male is estimated to be in the middle of this short period.
    age_tridactylophagus_tartari_days = 3
    age_tridactylophagus_tartari_unit = "days"

    # Species 2: Periplaneta americana (American cockroach)
    # Males reach sexual maturity and are capable of mating within 24 hours of eclosion.
    age_periplaneta_americana_days = 1
    age_periplaneta_americana_unit = "day"

    print("Estimated average age since eclosion for a mating male:")
    print(f"- Tridactylophagus tartari: {age_tridactylophagus_tartari_days} {age_tridactylophagus_tartari_unit}")
    print(f"- Periplaneta americana: {age_periplaneta_americana_days} {age_periplaneta_americana_unit}")

estimate_mating_age()