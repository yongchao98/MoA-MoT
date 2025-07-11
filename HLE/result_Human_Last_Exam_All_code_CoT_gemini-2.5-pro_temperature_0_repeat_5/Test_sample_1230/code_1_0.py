def analyze_lifestyle_demands():
    """
    This function models the demands on Charles to identify the area
    needing adjustment first.
    """
    # Assign difficulty scores to each component of Charles's plan.
    # The diet is broken down into its three restrictive components.
    limit_animal_protein_difficulty = 3
    limit_fat_difficulty = 3
    limit_carbs_difficulty = 4  # Limiting carbs often has the biggest initial impact on energy.
    
    new_exercise_difficulty = 3
    travel_difficulty = 2

    # Calculate the total difficulty for the diet
    total_diet_difficulty = limit_animal_protein_difficulty + limit_fat_difficulty + limit_carbs_difficulty
    
    # Calculate the overall total demand score
    total_demand = total_diet_difficulty + new_exercise_difficulty + travel_difficulty

    # Create a dictionary to hold the scores for comparison
    demands = {
        "Diet": total_diet_difficulty,
        "Exercise": new_exercise_difficulty,
        "Travel": travel_difficulty
    }

    # Find the component with the highest difficulty score
    primary_struggle = max(demands, key=demands.get)

    print("Analyzing the demands of the new lifestyle plan:")
    print(f"Total Demand Score = (Diet Difficulty) + (Exercise Difficulty) + (Travel Difficulty)")
    
    # The final equation with each number outputted
    print(f"The final equation is: {total_demand} = ({limit_animal_protein_difficulty} + {limit_fat_difficulty} + {limit_carbs_difficulty}) + {new_exercise_difficulty} + {travel_difficulty}")
    
    print(f"\nBreakdown of Demand Scores:")
    print(f"- Diet: {total_diet_difficulty}")
    print(f"- Exercise: {new_exercise_difficulty}")
    print(f"- Travel: {travel_difficulty}")

    print(f"\nConclusion:")
    print(f"The '{primary_struggle}' component contributes the most to the overall demand.")
    print("Therefore, Charles should first adjust the highly restrictive nature of his diet to make it more sustainable.")

analyze_lifestyle_demands()