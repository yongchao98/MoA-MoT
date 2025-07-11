def solve_riddle():
    """
    Solves the medical riddle by extracting numbers from the context
    and connecting them to the nutritional clue.
    """
    # 1. Extract the two numbers from the case history.
    patient_age = 29
    history_duration = 8

    # 2. The "bean salad" clue points to protein, which is made of amino acids.
    # 3. Perform the calculation.
    amino_acids_count = patient_age - history_duration

    # 4. Print the explanation.
    print("The importance of the new food, which tastes like 'bean salad', lies in its nutritional content.")
    print("Beans are an excellent source of protein, and proteins are made from amino acids.")
    print("The key to understanding its specific importance is hidden in the numbers from the patient's history.")
    print(f"By subtracting the duration of the history ({history_duration}) from the patient's age ({patient_age}), we get the number {amino_acids_count}.")
    print(f"This number, {amino_acids_count}, is the total number of proteinogenic amino acids required for life.")
    print("\nThe equation is:")
    
    # 5. Output the final equation with each number.
    print(f"{patient_age} - {history_duration} = {amino_acids_count}")

solve_riddle()