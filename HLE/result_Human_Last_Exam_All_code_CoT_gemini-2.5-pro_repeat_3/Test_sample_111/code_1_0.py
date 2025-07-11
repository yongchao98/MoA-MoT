def solve_historical_question():
    """
    Analyzes historical facts to determine Geoffrey Chaucer's location
    when Blanche of Lancaster died.
    """

    # Known facts
    blanche_death_date = "September 12, 1368"
    blanche_death_location = "England"
    chaucer_status_1368 = "Esquire in the English royal court"
    chaucer_activity_1368 = "Traveled to the continent on a diplomatic mission"

    # The problem: uncertainty in the timeline
    uncertainty = "The exact dates of Chaucer's 1368 diplomatic mission are not known."

    print("Investigating the location of Geoffrey Chaucer on the date of Blanche of Lancaster's death.")
    print(f"1. Blanche of Lancaster died on: {blanche_death_date} in {blanche_death_location}.")
    print(f"2. In the year 1368, Geoffrey Chaucer was an {chaucer_status_1368}.")
    print(f"3. During 1368, Chaucer's documented activity was: '{chaucer_activity_1368}'.")
    print(f"4. Point of uncertainty: {uncertainty}")
    print("\nConclusion:")
    print("Because Chaucer's precise return date from the continent is not recorded,")
    print("it is impossible to know for certain if he was in England or abroad when Blanche died.")
    print("Therefore, his exact location at that specific time is considered unknown by historians.")

solve_historical_question()