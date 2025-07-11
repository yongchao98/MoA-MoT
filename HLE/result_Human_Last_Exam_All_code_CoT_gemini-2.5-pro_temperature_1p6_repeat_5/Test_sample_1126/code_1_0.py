import numpy as np

def solve_fish_larvae_problem():
    """
    Analyzes fish larvae settlement data to determine the effect of CO2 levels
    and selects the best conclusion from a list of choices.
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

    print("Step 1 & 2: Calculating Average Time Spent Near Speaker\n")
    averages = {}
    for soundscape, conditions in data.items():
        averages[soundscape] = {}
        print(f"--- {soundscape} ---")
        for condition, percentages in conditions.items():
            avg = np.mean(percentages)
            averages[soundscape][condition] = avg
            print(f"Average for {condition}: {avg:.2f}%")
        print()

    print("Step 3: Analyzing Natural Habitat Preference (Control Condition)")
    print("A percentage > 50% indicates attraction, < 50% indicates avoidance.")
    print(f"- Tropical Estuarine (Control): {averages['Tropical estuarine soundscape']['control']:.2f}%. This clear attraction suggests it is a natural habitat.")
    print(f"- Temperate Reef (Control): {averages['Temperate reef soundscape']['control']:.2f}%. This is close to neutral, but could be a less preferred habitat.")
    print(f"- White Noise (Control): {averages['White noise']['control']:.2f}%. This suggests avoidance of an irrelevant/aversive sound.\n")

    print("Step 4: Analyzing the Effect of High CO2")
    print("Comparing 'Co2' averages to 'control' averages reveals behavioral changes.")
    print("- For the Tropical Estuarine soundscape, the preference under control conditions (56.33%) is lost and reverses to avoidance under high CO2 (43.67%). This is a significant disturbance.")
    print("- For White Noise, the avoidance under control conditions (46.17%) reverses to a strong attraction under high CO2 (53.83%). This shows a failure to distinguish irrelevant sounds.")
    print("Conclusion: High CO2 levels disrupt the larvae's sensory abilities, impairing their ability to use sound cues for settlement.\n")

    print("Step 5: Evaluating the Answer Choices")
    print("Based on the analysis:")
    print("- A is incorrect. Settlement behavior in the temperate reef is altered, not 'as efficiently'.")
    print("- B is incorrect. The data most strongly supports the tropical estuarine soundscape as the natural habitat, not just the temperate reef.")
    print("- D is too general. While true, it doesn't capture the specific findings of the experiment.")
    print("- E is incorrect. High CO2 clearly disturbs settlement behavior related to the tropical estuarine sound.")
    print("- F is incorrect. It rightly states the disturbance for the tropical estuarine but wrongly claims the temperate reef is the sole natural habitat.")
    print("- H is plausible, but C is more specific about the nature of the disturbance.")
    print("\n- C is the best answer. It correctly identifies that settlement in the tropical estuarine habitat will be inefficient under high CO2. It also reasonably assumes that both soundscapes could be considered natural habitats, which is consistent with the experimental design.")
    
    print("\nFinal Conclusion:")
    print("The data shows that under normal conditions, larvae are attracted to the Tropical Estuarine soundscape (avg = 56.33%). Under high CO2, this attraction is lost (avg = 43.67%). This means the fish will not settle in the tropical estuarine habitat as efficiently.")
    print("Therefore, choice C is the most accurate conclusion.")
    print("C. The natural habitats of the examined fish are the tropical estuarine and the temperate reef. At the Co2 level predicted for the year 2100, the fish will not settle in the tropical estuarine efficiently as in the year 2024.")

    final_answer = 'C'
    print(f"\n<<<{final_answer}>>>")

solve_fish_larvae_problem()