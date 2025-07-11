import sys

def solve_historical_question():
    """
    This function determines Geoffrey Chaucer's location
    at the time of Blanche of Lancaster's death.
    """
    # Step 1: Establish the date and location of Blanche of Lancaster's death.
    year_of_death = 1368
    location_of_death = "England (at Bolingbroke Castle, Lincolnshire)"

    print(f"Fact 1: Blanche of Lancaster died in the year {year_of_death} in {location_of_death}.")

    # Step 2: Establish Geoffrey Chaucer's status and location at that time.
    chaucer_role_in_1367 = "made an esquire in the household of King Edward III"
    royal_court_location = "England"
    first_major_poem = "The Book of the Duchess"
    poem_subject = "an elegy for Blanche of Lancaster"

    print(f"Fact 2: By {year_of_death}, Chaucer was serving as {chaucer_role_in_1367}.")
    print(f"Fact 3: The royal court he served was based in {royal_court_location}.")
    print(f"Fact 4: Chaucer's first major poem, '{first_major_poem}', was written in her memory, suggesting a close connection to the court at the time.")
    print("Fact 5: His famous diplomatic journeys to Italy and France happened in the 1370s, after her death.")

    # Step 3: Conclude Chaucer's location.
    conclusion = "Based on these facts, Geoffrey Chaucer was in England when Blanche of Lancaster died."
    print("\nConclusion: " + conclusion)
    
    # Step 4: Identify the corresponding answer choice.
    answer_choice = 'C'
    print(f"\nThe correct answer choice is '{answer_choice}'.")

solve_historical_question()