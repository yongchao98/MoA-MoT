def solve_insect_riddle():
    """
    Solves the riddle about the feeding habits of an adult male Micromalthidae.
    """
    # Step 1: Store biological facts about Micromalthidae feeding habits.
    feeding_habits = {
        'larva': 'Fungus-infested, decaying wood. Some specialized larvae also exhibit matriphagy (eating their mother).',
        'adult_female': 'Believed to be non-feeding, with a very short lifespan.',
        'adult_male': 'Nothing. Adult males are extremely rare, short-lived, have non-functional (vestigial) mouthparts, and do not feed.'
    }

    # Step 2: Identify the subject of the question.
    subject_life_stage = 'adult_male'

    # Step 3: Retrieve the relevant biological fact.
    fact = feeding_habits.get(subject_life_stage, 'Unknown life stage.')

    # Step 4: Present the reasoning and the answer.
    print("Analyzing the feeding habits of Micromalthidae:")
    print("-" * 50)
    print(f"The question concerns the diet of an: adult male.")
    print(f"The established biological fact for this life stage is:")
    print(f"    - {fact}")
    print("-" * 50)
    print("Comparing this fact to the answer choices:")
    print("A. Its mother - This is a larval behavior (matriphagy), not adult male behavior.")
    print("B. Fungus - This is food for larvae.")
    print("C. Decaying wood - This is food for larvae.")
    print("D. Cellulose - A component of wood, which is food for larvae.")
    print("E. Nothing - This matches the fact that adult males are non-feeding.")
    
    print("\nConclusion: Upon its death, the only thing this individual will have fed on is nothing, as adult males do not eat.")
    print("The correct answer is E.")

solve_insect_riddle()