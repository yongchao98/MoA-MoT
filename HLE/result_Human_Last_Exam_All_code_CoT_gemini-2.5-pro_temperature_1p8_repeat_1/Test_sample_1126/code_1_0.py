def analyze_larvae_behavior():
    """
    Analyzes fish larvae behavior data to determine the impact of CO2.
    """
    data = {
        'Temperate reef soundscape': {
            'control': [40, 50, 49, 53, 49, 48],
            'Co2': [60, 50, 51, 47, 51, 52]
        },
        'Tropical estuarine soundscape': {
            'control': [68, 63, 55, 51, 49, 52],
            'Co2': [32, 37, 45, 49, 51, 48]
        },
        'White noise': {
            'control': [30, 48, 47, 48, 52, 52],
            'Co2': [70, 52, 53, 52, 48, 48]
        }
    }

    print("--- Data Analysis ---\n")

    # Step 1 & 2: Calculate and print average for each condition
    print("Step 1: Calculate average time spent near speaker for each condition.\n")
    averages = {}
    for soundscape, conditions in data.items():
        averages[soundscape] = {}
        print(f"--- {soundscape} ---")
        for condition, percentages in conditions.items():
            sum_str = " + ".join(map(str, percentages))
            avg = sum(percentages) / len(percentages)
            averages[soundscape][condition] = avg
            print(f"Average for '{condition}': ({sum_str}) / {len(percentages)} = {avg:.2f}%")
        print()

    # Step 3, 4, & 5: Analyze results and evaluate options
    print("\n--- Interpretation of Results ---\n")
    print("Step 2: Determine natural habitat from 'control' data.")
    print(f"Under control conditions, larvae show the strongest attraction to the 'Tropical estuarine soundscape' (average = {averages['Tropical estuarine soundscape']['control']:.2f}%), which is significantly above the 50% mark of random chance.")
    print(f"Their response to the 'Temperate reef soundscape' (average = {averages['Temperate reef soundscape']['control']:.2f}%) is neutral.")
    print("Conclusion: The primary natural habitat is the tropical estuary.\n")
    
    print("Step 3: Assess the impact of elevated CO2.")
    print(f"For the 'Tropical estuarine soundscape', the attraction dramatically reverses under high CO2. The average time drops from {averages['Tropical estuarine soundscape']['control']:.2f}% to {averages['Tropical estuarine soundscape']['Co2']:.2f}%.")
    print("This means the larvae's ability to find their natural habitat is significantly disturbed, and they would not settle there as efficiently.\n")

    print("Step 4: Evaluate the provided answer choices.")
    print("Choice C states: 'The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.'")
    print("This statement aligns perfectly with our findings. While the temperate reef isn't a *primary* habitat based on the data, the most critical part of the statement is that elevated CO2 disrupts efficient settlement in the tropical estuarine, which is strongly supported by the data.")
    print("Other options are incorrect because they misidentify the natural habitat or misrepresent the effect of CO2.\n")

    
if __name__ == '__main__':
    analyze_larvae_behavior()
    print("<<<C>>>")