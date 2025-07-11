def solve_opera_riddle():
    """
    This function solves a multi-step opera trivia question.
    """

    # 1. Define the key years and numbers from the clues.
    caruso_last_huguenots_year = 1915
    time_gap = 70
    nyc_revival_year = 1991

    # 2. Perform the calculation to establish the required timeframe for the revival.
    # The revival had to occur after this calculated year.
    min_revival_year = caruso_last_huguenots_year + time_gap
    
    # Calculate the actual gap in years to confirm it fits the clue.
    actual_time_passed = nyc_revival_year - caruso_last_huguenots_year

    # 3. Identify the key details of the opera and the performance.
    opera_title = "Les Huguenots"
    bass_role_name = "Marcel"
    bass_singer_name = "Nicolai Ghiaurov"

    # 4. Print the reasoning and the solution step-by-step.
    print(f"The opera is identified as '{opera_title}'.")
    print(f"Enrico Caruso's last performance in this opera at the Met was in {caruso_last_huguenots_year}.")
    print("The revival had to be more than 70 years later.")
    
    print("\n--- Timeline Calculation ---")
    print(f"The last performance year was {caruso_last_huguenots_year}.")
    print(f"The required time gap was {time_gap} years.")
    print(f"Equation: {caruso_last_huguenots_year} + {time_gap} = {min_revival_year}")
    print(f"Therefore, the production had to be after {min_revival_year}.")
    
    print(f"\nA major NYC revival occurred in {nyc_revival_year}, which is {actual_time_passed} years later, fitting the criteria.")
    
    print("\n--- Final Answer ---")
    print(f"In the {nyc_revival_year} New York production of '{opera_title}', the bass role of {bass_role_name} was sung by:")
    print(bass_singer_name)

solve_opera_riddle()