def analyze_lifestyle_change(components):
    """
    Analyzes components of a lifestyle change to identify the area
    with the highest demand, which is the most likely source of struggle.

    Args:
        components (dict): A dictionary where keys are the component names
                           and values are their associated demand scores.
    
    Returns:
        str: The name of the component that needs adjustment first.
    """
    if not components:
        return "No components to analyze."

    # Find the component with the highest demand score
    primary_struggle_area = max(components, key=components.get)
    return primary_struggle_area

# Assign a "demand score" to each part of Charles's plan. A higher score
# indicates a more difficult component to maintain.
# The diet is extremely restrictive, so it gets the highest score.
charles_plan_components = {
    "Highly Restrictive Diet (limiting protein, fat, and carbs)": 10,
    "Demanding Lifestyle (frequent travel)": 5,
    "Planned New Exercise Routine": 4
}

# Get the component that is the primary source of struggle
first_adjustment_needed = analyze_lifestyle_change(charles_plan_components)
highest_score = charles_plan_components[first_adjustment_needed]

# --- Output the analysis ---
print("Step-by-Step Analysis of Charles's Lifestyle Demands:")
print("=====================================================")

# Print each component and its score to show the 'equation'
for component, score in charles_plan_components.items():
    print(f"Demand Score for '{component}': {score}")

print("\nConclusion:")
print(f"The component with the highest demand score is '{first_adjustment_needed}' with a score of {highest_score}.")
print("This indicates it is the most significant source of his struggle.")

print("\nRecommendation:")
print("Charles should first adjust his diet. A diet that severely limits protein, fat, and carbohydrates all at once is extremely difficult to sustain, especially with a demanding travel schedule. He needs to create a more balanced and less restrictive eating plan before incorporating a new exercise routine. This will create a sustainable foundation for long-term success.")