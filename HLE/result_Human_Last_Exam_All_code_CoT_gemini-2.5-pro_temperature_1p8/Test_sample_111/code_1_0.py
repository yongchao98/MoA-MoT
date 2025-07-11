def find_chaucer_location():
    """
    This function determines Geoffrey Chaucer's location
    at the time of Blanche of Lancaster's death by analyzing historical data.
    """

    # Step 1: Establish the year of Blanche of Lancaster's death.
    year_of_death = 1368

    print(f"Fact: Blanche of Lancaster died in the year {year_of_death}, likely from the plague.")

    # Step 2: Access historical data about Chaucer's activities around that year.
    # By 1367, Chaucer was a member of the English royal court as an esquire to King Edward III.
    # His famous poem, 'The Book of the Duchess', is an elegy for Blanche of Lancaster,
    # and it is believed to have been written between 1368 and 1372.
    # His first documented diplomatic mission to Italy occurred later, in 1372.
    
    chaucer_status_1367 = "Esquire in the royal court of King Edward III in England."
    chaucer_first_italian_journey = 1372

    print(f"Fact: In 1367, Geoffrey Chaucer's position was: '{chaucer_status_1367}'.")
    print(f"Fact: Chaucer's first recorded trip to Italy was not until {chaucer_first_italian_journey}.")

    # Step 3: Conclude based on the evidence.
    conclusion = (
        "Based on these facts, Chaucer was serving the royal court in England "
        f"in {year_of_death}. His poem mourning Blanche's death further "
        "solidifies his presence within the English court circle at that time."
    )
    
    print("\nConclusion: " + conclusion)
    print("\nAnswer: C. England")

find_chaucer_location()