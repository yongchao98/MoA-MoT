import textwrap

def analyze_tree_ring_data():
    """
    Analyzes the factors influencing the 13C isotope ratio in tree rings.
    The core question is about a *declining* 13C ratio from 1886-1990.
    """

    question = "Which of the factors below was the predominant factor influencing the declining 13C ratio observed in the tree ring data over this period?"
    
    options = {
        'A': "An increase in tree ring thickness as the tree matures",
        'B': "Periods of drought",
        'C': "Increased photosynthetic reserves of starch fueling tree growth",
        'D': "Thinning earlywood tree ring proportion",
        'E': "Changes in the SE Asia monsoon"
    }

    print("--- Analyzing the Predominant Factor for Declining Tree Ring 13C Ratio (1886-1990) ---")
    print("\n[Background Information]")
    print(textwrap.fill(
        "A decline in the 13C/12C ratio (δ13C) in a tree ring can be caused by two primary phenomena: "
        "1) A decline in the δ13C of the source atmospheric CO2. "
        "2) An increase in the tree's discrimination against 13C during photosynthesis, often due to environmental factors like increased water availability.", 80
    ))
    
    print("\n[Evaluation of Answer Choices]\n")

    # Analysis of Option B (Drought)
    print("Analysis of Option B: Periods of drought")
    print(textwrap.fill(
        "Drought stress causes a tree's stomata (leaf pores) to close. This reduces the CO2 supply inside the leaf, forcing the photosynthetic enzymes to be less selective. This leads to an INCREASE in the 13C ratio, which is the opposite of the observed decline. Therefore, B is incorrect.", 80
    ))
    print("-" * 80)
    
    # Analysis of Options A, C, D (Local/Physiological Factors)
    print("Analysis of Options A, C, and D: Tree-specific factors")
    print(textwrap.fill(
        "Factors like tree ring thickness (A), starch reserves (C), and earlywood proportion (D) are related to the individual tree's age, health, and internal physiology. While they can cause minor variations, they are not considered large-scale 'predominant' drivers for a consistent, century-long trend across a region.", 80
    ))
    print("-" * 80)

    # Analysis of Option E (SE Asia Monsoon)
    print("Analysis of Option E: Changes in the SE Asia monsoon")
    print(textwrap.fill(
        "The SE Asia monsoon is a massive climate system that controls regional water availability. A long-term trend towards a stronger or wetter monsoon season means trees experience less water stress. With ample water, their stomata can remain open wider, allowing for greater discrimination against the heavier 13C isotope. This leads to a DECLINING 13C ratio in the tree rings. Of the choices provided, a large-scale climate driver like the monsoon is the most plausible 'predominant factor' capable of producing such a long-term, directional change.", 80
    ))
    print("-" * 80)

    # Conclusion
    print("\n[Conclusion]")
    print(textwrap.fill(
        "While the global decline in atmospheric δ13C due to fossil fuel burning (the Suess Effect) is the ultimate underlying cause, it is not an option. Among the choices given, Changes in the SE Asia monsoon (E) is the best fit as a predominant environmental factor that can cause a long-term decline in the 13C ratio.", 80
    ))
    
    final_choice = 'E'
    print(f"\nThe best answer from the list is: {final_choice}")

# Run the analysis
analyze_tree_ring_data()