import textwrap

def analyze_tree_ring_factors():
    """
    This script evaluates the potential factors influencing the 13C isotope ratio
    in Chinese pine trees between 1886-1990, as described in the problem.
    The goal is to identify the predominant factor for the observed decline.
    """

    # --- Scientific Principle ---
    # A declining 13C ratio in tree rings indicates that the tree was able to
    # discriminate more against the heavier 13C isotope during photosynthesis.
    # This typically happens under conditions of low water stress (i.e., wetter, more humid conditions),
    # which allows the tree's stomata (pores) to remain open.

    print("Analyzing the predominant factor for a declining 13C ratio in tree rings.\n")

    factors = {
        'A': "An increase in tree ring thickness as the tree matures.",
        'B': "Periods of drought.",
        'C': "Increased photosynthetic reserves of starch fueling tree growth.",
        'D': "Thinning earlywood tree ring proportion.",
        'E': "Changes in the SE Asia monsoon."
    }

    # --- Evaluation Logic ---
    # We will loop through each factor and assess its likely impact on the 13C ratio.

    best_factor = None
    print("--- Factor Evaluation ---")

    # Analysis for Factor B
    impact_b = "INCREASES 13C ratio (contradicts observation)."
    print(f"Factor B: {factors['B']}")
    print(f"Analysis: Drought causes water stress, forcing stomata to close. This reduces discrimination against 13C, causing the ratio to INCREASE.")
    print(f"Result: Incompatible with a declining ratio.\n")

    # Analysis for Factor E
    impact_e = "DECREASES 13C ratio (matches observation)."
    print(f"Factor E: {factors['E']}")
    print(textwrap.fill(
        f"Analysis: The monsoon is the primary driver of moisture in Southeast Asia. A strengthening monsoon over time leads to wetter conditions and lower water stress. This allows stomata to open, increasing discrimination against 13C and causing the ratio to DECREASE.",
        width=80
    ))
    print("Result: Plausible primary driver for a long-term decline.\n")

    # Briefly address other factors
    print("Brief analysis of other factors (A, C, D):")
    print(textwrap.fill(
        "Factors A, C, and D describe physiological or secondary effects within the tree. They are not considered primary environmental drivers capable of causing a consistent, predominant trend over a 100-year period across a region when compared to a major climatic system like the monsoon.",
        width=80
    ))
    
    best_factor = 'E'

    print("\n--- Conclusion ---")
    print(f"The most likely predominant factor is '{factors[best_factor]}'.")
    print("It is the only option that represents a large-scale climatic driver whose change would directly cause a long-term decline in the 13C ratio.")

analyze_tree_ring_factors()
<<<E>>>