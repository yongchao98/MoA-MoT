def find_chaucer_location():
    """
    This function uses historical data to determine Chaucer's likely location
    at the time of Blanche of Lancaster's death.
    """

    # Historical Facts
    death_of_blanche_year = 1368
    chaucer_role = "Esquire in the English royal court"
    record_of_travel_1368 = "No surviving records show Chaucer on a diplomatic mission abroad in the autumn of 1368."
    commemorative_work = "The Book of the Duchess, written in England for an English audience shortly after her death."

    # Step-by-step reasoning
    print("Step 1: Establish the date of Blanche of Lancaster's death.")
    print(f"Blanche of Lancaster died in the year {death_of_blanche_year}.")
    print("\nStep 2: Review Geoffrey Chaucer's status at that time.")
    print(f"In {death_of_blanche_year}, Chaucer's primary role was as an '{chaucer_role}'.")
    print(f"His home and place of work were in England.")
    print(f"Fact Check: {record_of_travel_1368}")

    print("\nStep 3: Draw a conclusion based on the evidence.")
    print("Given that Chaucer was in the service of the English court and not documented to be abroad, it is historically accepted that he was in England.")
    print(f"His poem, '{commemorative_work}', further supports this conclusion.")

    # Final Answer
    final_location = "England"
    print(f"\nConclusion: The evidence strongly indicates Geoffrey Chaucer was in {final_location}.")


# Run the analysis
find_chaucer_location()