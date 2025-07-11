def find_chaucers_location():
    """
    This function uses historical data to determine Chaucer's location
    at the time of Blanche of Lancaster's death.
    """
    # Fact 1: Blanche of Lancaster's death
    blanche_death_year = 1368
    blanche_death_location = "England"  # She died at Tutbury Castle, Staffordshire

    # Fact 2: Geoffrey Chaucer's status in that year
    chaucer_primary_location_1368 = "England"
    chaucer_role_1368 = "Esquire in the royal household of King Edward III"

    # Logical Deduction:
    # If Chaucer's primary base of operations was England and Blanche died in England,
    # it is overwhelmingly likely that he was in England at the time.
    # His famous poem, "The Book of the Duchess," was an elegy for her,
    # further cementing his connection to the English court's mourning process.

    print(f"Year of Blanche of Lancaster's death: {blanche_death_year}")
    print(f"Location of Blanche of Lancaster's death: {blanche_death_location}")
    print(f"Chaucer's primary location in {blanche_death_year}: {chaucer_primary_location_1368}")
    
    if chaucer_primary_location_1368 == blanche_death_location:
        conclusion = "England"
        print(f"\nConclusion: Based on the historical record, Geoffrey Chaucer was in {conclusion} when Blanche of Lancaster died.")
    else:
        # This path is not historically supported, but included for completeness.
        conclusion = "Unknown"
        print(f"\nConclusion: Location cannot be definitively determined, but evidence points away from {chaucer_primary_location_1368}.")

find_chaucers_location()