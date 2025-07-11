def analyze_fish_behavior():
    """
    Analyzes fish larvae behavior data to determine the effect of CO2.
    The function calculates the average time spent near speakers for each condition
    and then evaluates the given hypotheses.
    """
    data = {
        "Temperate reef soundscape": {
            "control": [40, 50, 49, 53, 49, 48],
            "co2": [60, 50, 51, 47, 51, 52]
        },
        "Tropical estuarine soundscape": {
            "control": [68, 63, 55, 51, 49, 52],
            "co2": [32, 37, 45, 49, 51, 48]
        },
        "White noise": {
            "control": [30, 48, 47, 48, 52, 52],
            "co2": [70, 52, 53, 52, 48, 48]
        }
    }

    print("--- Data Analysis ---\n")

    # Calculate and print averages for each group
    for soundscape, conditions in data.items():
        print(f"Analyzing: {soundscape}")
        for condition, percentages in conditions.items():
            num_days = len(percentages)
            total = sum(percentages)
            average = total / num_days
            # Create a string of the numbers being added for the print statement
            equation_str = " + ".join(map(str, percentages))
            print(f"  - {condition.upper()} Average: ({equation_str}) / {num_days} = {average:.2f}%")
        print("-" * 20)

    print("\n--- Conclusion ---\n")
    print("1. Identify Natural Habitat (from Control data):")
    print("   - Fish show strong attraction to the 'Tropical estuarine soundscape' (avg > 50%), indicating it's their preferred natural habitat.")
    print("   - They show indifference or weak attraction to the 'Temperate reef soundscape' (avg ≈ 48%).")
    print("   - They show repulsion or indifference to 'White noise' (avg < 50%).")
    
    print("\n2. Analyze Effect of High CO2:")
    print("   - For the 'Tropical estuarine soundscape', attraction drops from 56.33% to 43.67%. The larvae no longer find their habitat attractive.")
    print("   - This means at the CO2 level for 2100, the fish will not settle in the tropical estuarine as efficiently as they do now.")

    print("\n3. Evaluate Answer Choices:")
    print("   - A & B & E & F are incorrect because they misidentify the natural habitat or wrongly describe the effects.")
    print("   - H is too general. D is true but not the main point.")
    print("   - C correctly identifies that settlement in the tropical estuarine (the preferred habitat) will not be efficient under high CO2, which is the most significant finding.")
    print("   - Choice C assumes both are potential natural habitats, and focuses on the key change: the loss of attraction to the primary habitat.")

analyze_fish_behavior()
<<<C>>>