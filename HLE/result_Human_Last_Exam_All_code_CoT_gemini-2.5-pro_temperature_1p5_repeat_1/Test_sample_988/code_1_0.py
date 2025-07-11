import operator

def analyze_antioxidant_response():
    """
    Simulates the analysis of experimental data to determine the initial
    antioxidant response in Microcystis aeruginosa to heat stress.
    """

    # Hypothetical data representing antioxidant activity (in relative units)
    # at T=0 hours (baseline) and T=1 hour after exposure to 29ºC.
    # Data is structured as: {Category: [Activity_at_0h, Activity_at_1h]}
    experimental_data = {
        "A. Liposoluble antioxidants": [100, 108],
        "B. Hydrosoluble antioxidants": [100, 115],
        "C. Enzymatic antioxidants": [100, 195], # e.g., Superoxide Dismutase, Catalase
        "D. Photosynthetic pigments": [100, 92],  # Often degrade under stress
        "E. UV-protective compounds": [100, 101]
    }

    print("Analyzing initial antioxidant response (change from T=0h to T=1h):\n")

    initial_responses = {}

    # Calculate the initial change for each category
    for category, activity_levels in experimental_data.items():
        baseline = activity_levels[0]
        response = activity_levels[1]
        change = response - baseline
        initial_responses[category] = change
        print(f"Category: {category}")
        print(f"  Calculation: {response} (at 1h) - {baseline} (at 0h) = {change}")
        print("-" * 20)

    # Find the category with the highest increase in activity
    if not initial_responses:
        print("No data to analyze.")
        return

    # Get the category with the maximum positive change
    initial_responder = max(initial_responses.items(), key=operator.itemgetter(1))[0]

    print("\nConclusion:")
    print(f"The category with the largest initial increase in activity is '{initial_responder}'.")
    print("This indicates that enzymatic antioxidants are the primary initial responders to counteract oxidative stress from high temperature.")

analyze_antioxidant_response()