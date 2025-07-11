import textwrap

def analyze_tree_ring_data():
    """
    This script analyzes the provided multiple-choice question about tree ring isotopes
    and determines the most plausible answer based on scientific principles.
    """

    # Key numbers from the problem description
    start_year = 1886
    end_year = 1990
    isotope_mass_number = 13

    question = f"A study of Chinese pine tree rings from {start_year}-{end_year} AD showed a declining ratio of the {isotope_mass_number}C isotope. The task is to identify the predominant factor for this decline from the given options."

    print("Step 1: Understanding the core observation.")
    print(textwrap.fill("The key observation is a 'declining 13C ratio'. In isotope science, this means the sample is becoming more depleted in the heavy 13C isotope relative to the lighter 12C isotope. This can be caused by changes in the source CO2 or by plant physiological responses to environmental conditions like water availability.", 80))
    print("-" * 80)

    options = {
        'A': "An increase in tree ring thickness as the tree matures",
        'B': "Periods of drought",
        'C': "Increased photosynthetic reserves of starch fueling tree growth",
        'D': "Thinning earlywood tree ring proportion",
        'E': "Changes in the SE Asia monsoon"
    }

    rationales = {
        'A': "Incorrect. This describes an age-related trend. While it can influence the record, it is not the predominant environmental factor explaining a consistent trend over a century.",
        'B': "Incorrect. Drought stress causes a plant's stomata to close. This leads to *less* discrimination against the heavy 13C isotope, resulting in an *increase* in the 13C ratio, which is the opposite of the observed trend.",
        'C': "Incorrect. Starch reserves reflect recent photosynthesis. This internal allocation strategy is not the primary driver of a multi-decadal environmental trend.",
        'D': "Incorrect. The proportion of earlywood to latewood is typically a *result* of climatic conditions, not the ultimate cause itself.",
        'E': "Correct. The SE Asia monsoon is the main driver of regional climate, especially precipitation. A long-term strengthening of the monsoon would increase water availability. With more water, trees can keep their stomata open wider and for longer, allowing for greater enzymatic discrimination against the heavier 13C isotope. This results in a declining 13C ratio in the wood formed each year."
    }

    print("Step 2: Evaluating each potential factor.")
    for key in options:
        print(f"\n[Option {key}] {options[key]}")
        print(textwrap.fill(f"Analysis: {rationales[key]}", 80, initial_indent="  ", subsequent_indent="  "))

    correct_answer = 'E'

    print("-" * 80)
    print("Step 3: Final Conclusion.")
    # Although there's no mathematical equation to solve, we can display the key numbers
    # involved in the problem's context as requested.
    print(f"The period under consideration is from the year {start_year} to {end_year}.")
    print(f"The analysis concerns the carbon isotope with mass number {isotope_mass_number}.")
    print(f"\nBased on the analysis, the predominant factor is Option {correct_answer}: '{options[correct_answer]}'.")

analyze_tree_ring_data()

print("\n<<<E>>>")