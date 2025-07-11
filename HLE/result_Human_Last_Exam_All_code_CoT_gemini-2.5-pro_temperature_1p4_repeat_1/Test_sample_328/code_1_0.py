def estimate_mating_age():
    """
    This function provides the best average estimate of the age since eclosion
    for a male Tridactylophagus tartari and a male Periplaneta americana to be observed mating.
    """

    # Species 1: Tridactylophagus tartari (a parasitic digger wasp)
    # Research indicates that digger wasps often have a maturation feeding period
    # after eclosion that lasts a few days before they are ready to mate.
    wasp_species = "Tridactylophagus tartari"
    wasp_age_estimate_days = 3
    wasp_age_unit = "days"

    # Species 2: Periplaneta americana (American cockroach)
    # Research indicates that male American cockroaches can become sexually
    # mature and mate very quickly after their final molt, with mating
    # possible as early as the first day.
    cockroach_species = "Periplaneta americana"
    cockroach_age_estimate_days = 1
    cockroach_age_unit = "day"

    print(f"Based on entomological data:")
    print(f"The best average age estimate for a male {wasp_species} is {wasp_age_estimate_days} {wasp_age_unit} since eclosion.")
    print(f"The best average age estimate for a male {cockroach_species} is {cockroach_age_estimate_days} {cockroach_age_unit} since eclosion.")

estimate_mating_age()
<<<D>>>