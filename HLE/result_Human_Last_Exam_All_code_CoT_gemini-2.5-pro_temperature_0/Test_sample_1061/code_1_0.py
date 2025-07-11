import textwrap

def analyze_tree_ring_data():
    """
    Analyzes the potential factors for the declining 13C ratio in tree rings
    from 1886-1990 and prints a step-by-step evaluation.
    """
    observation = "A declining 13C ratio in Chinese pine tree rings over the period of 1886-1990 AD."
    options = {
        'A': "An increase in tree ring thickness as the tree matures",
        'B': "Periods of drought",
        'C': "Increased photosynthetic reserves of starch fueling tree growth",
        'D': "Thinning earlywood tree ring proportion",
        'E': "Changes in the SE Asia monsoon"
    }

    print(f"Observation: {observation}\n")
    print("Evaluating the potential predominant factors:\n" + "="*40)

    # --- Evaluation Logic ---

    # Option A
    analysis_a = "This is a known 'juvenile effect'. While it can influence the 13C ratio, it is typically a factor in the early decades of growth and is unlikely to be the predominant driver of a continuous, 100+ year trend."
    print(f"1. Evaluating Option A: {options['A']}")
    print(textwrap.fill(f"   Analysis: {analysis_a}\n", width=80))

    # Option B
    analysis_b = "Drought causes water stress, forcing the tree to close its stomata. This limits the CO2 supply and reduces the tree's ability to discriminate against 13C, leading to an INCREASE in the 13C ratio. This is the opposite of the observed trend."
    print(f"2. Evaluating Option B: {options['B']}")
    print(textwrap.fill(f"   Analysis: {analysis_b}\n", width=80))

    # Option C
    analysis_c = "The use of stored starch is a secondary physiological process. The isotopic signature of the starch itself was set by primary factors (like climate and atmospheric CO2) when it was created. This is not a primary driver of a long-term trend."
    print(f"3. Evaluating Option C: {options['C']}")
    print(textwrap.fill(f"   Analysis: {analysis_c}\n", width=80))

    # Option D
    analysis_d = "Changes in the proportion of earlywood to latewood are typically a response to the same environmental conditions that influence isotope ratios, rather than being the primary cause of the isotopic trend itself."
    print(f"4. Evaluating Option D: {options['D']}")
    print(textwrap.fill(f"   Analysis: {analysis_d}\n", width=80))

    # Option E
    analysis_e = "The SE Asia monsoon is a major driver of regional climate, especially moisture availability. A long-term strengthening of the monsoon would lead to less water stress. This allows the tree's stomata to remain open, enhancing discrimination against 13C and causing a DECLINING 13C ratio. This aligns with the observation."
    print(f"5. Evaluating Option E: {options['E']}")
    print(textwrap.fill(f"   Analysis: {analysis_e}\n", width=80))

    print("="*40 + "\nConclusion:\n")
    
    conclusion_text = ("While the global 'Suess Effect' (declining atmospheric 13C from fossil fuels) is the ultimate global driver for this period, it is not an option. Among the choices provided, a large-scale climatic shift is the most plausible predominant factor. A stronger SE Asia monsoon leads to increased moisture and a lower 13C ratio.")
    print(textwrap.fill(conclusion_text, width=80))

    # The "final equation" is a representation of the contributing factors.
    print("\nFinal Equation of Influences:")
    print("Observed 13C Decline = Effect(Changes in SE Asia monsoon) + Effect(Global Atmospheric Change)")
    print("\nBased on the provided choices, the most significant factor is E.")


if __name__ == '__main__':
    analyze_tree_ring_data()
    print("\n<<<E>>>")