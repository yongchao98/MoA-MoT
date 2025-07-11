def solve_historical_question():
    """
    This script determines where Geoffrey Chaucer was when Blanche of Lancaster died
    by using established historical facts.
    """
    
    # Historical facts stored as variables
    year_of_death = 1368
    person_who_died = "Blanche of Lancaster"
    person_in_question = "Geoffrey Chaucer"
    
    # Chaucer's known location and status around the time of death
    chaucers_location = "England"
    chaucers_role = "an esquire in the royal household of King Edward III"
    commemorative_work = "The Book of the Duchess"
    
    # Print the logical deduction step-by-step
    print(f"Step 1: Determine the year of death for {person_who_died}.")
    # The 'final equation' will be represented by showing the key number.
    print(f"The year of death is = {year_of_death}")
    print("-" * 30)
    
    print(f"Step 2: Determine the location of {person_in_question} in the year {year_of_death}.")
    print(f"Historical records show that in {year_of_death}, Chaucer was serving as {chaucers_role}.")
    print(f"This service placed him firmly in {chaucers_location}.")
    print("-" * 30)

    print(f"Step 3: Corroborate with related works.")
    print(f"Chaucer's first major work, '{commemorative_work}', was an elegy written for Blanche of Lancaster, strongly tying him to the English court at the time of her death.")
    print("-" * 30)
    
    print("\nConclusion:")
    print(f"Based on the evidence, {person_in_question} was in {chaucers_location} when {person_who_died} died in {year_of_death}.")

solve_historical_question()