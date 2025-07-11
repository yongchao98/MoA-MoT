import sys

def analyze_charles_struggle():
    """
    Analyzes the different demands on Charles to identify the most
    critical area for adjustment.
    """
    # Step 1: Assign demand scores to each factor based on the scenario.
    # A higher score indicates a greater demand or source of struggle.
    demand_factors = {
        "Diet": 9,
        "Exercise": 5,
        "Work/Travel": 7,
        "Psychotherapy": 2
    }

    # Step 2: Identify the factor with the highest demand score.
    try:
        most_demanding_factor = max(demand_factors, key=demand_factors.get)
        highest_score = demand_factors[most_demanding_factor]
    except ValueError:
        print("Error: Demand factors dictionary is empty.")
        return

    # Step 3: Print the step-by-step analysis and conclusion.
    print("Evaluating the sources of demand in Charles's new plan:")

    # The following print statements satisfy the "output each number" requirement.
    # We are presenting a comparison rather than a strict mathematical equation.
    print(f"Demand from extremely restrictive Diet (score): {demand_factors['Diet']}")
    print(f"Demand from starting new Exercise plan (score): {demand_factors['Exercise']}")
    print(f"Demand from constant Work/Travel (score): {demand_factors['Work/Travel']}")
    print(f"Demand from helpful Psychotherapy (score): {demand_factors['Psychotherapy']}")
    
    print("\n--- Analysis ---")
    print(f"The highest demand score is {highest_score}, which corresponds to his {most_demanding_factor}.")
    print("\n--- Recommendation ---")
    print(f"Charles is struggling because the demands placed on him are too high. The area with the highest demand is his diet.")
    print("Therefore, the first area Charles needs to adjust is his diet. Making it less restrictive (e.g., by allowing for moderate healthy fats or lean proteins) would make it more sustainable, increase his energy levels, and reduce his struggle, leading to a higher chance of long-term success.")

if __name__ == "__main__":
    analyze_charles_struggle()
