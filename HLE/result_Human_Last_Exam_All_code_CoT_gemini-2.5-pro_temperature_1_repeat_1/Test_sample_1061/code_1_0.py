import textwrap

def analyze_tree_ring_data():
    """
    Analyzes the factors influencing the 13C ratio in tree rings to determine the correct answer.
    """

    question = "Which of the factors below was the predominant factor influencing the declining 13C ratio observed in the tree ring data over this period (1886-1990AD)?"
    
    options = {
        'A': "An increase in tree ring thickness as the tree matures",
        'B': "Periods of drought",
        'C': "Increased photosynthetic reserves of starch fueling tree growth",
        'D': "Thinning earlywood tree ring proportion",
        'E': "Changes in the SE Asia monsoon"
    }

    print("Analyzing the question: " + question)
    print("-" * 70)

    # Scientific background
    print("Scientific Background:")
    background = """
    A declining 13C ratio (δ13C) in tree rings primarily reflects two things:
    1. The isotopic composition of the source CO2 in the atmosphere. The 'Suess Effect' from burning 13C-depleted fossil fuels has caused a global decline in atmospheric δ13C since the industrial revolution.
    2. The plant's physiological response to its environment. Higher water availability allows the tree to open its stomata wider, increasing its discrimination against the heavier 13C isotope and resulting in a lower δ13C in its wood.
    """
    print(textwrap.dedent(background))
    print("-" * 70)

    # Analysis of each option
    print("Evaluation of Answer Choices:")
    
    # Option A
    analysis_a = "Incorrect. While age can have some effect, it is not the predominant driver of a large-scale, 100-year isotopic trend compared to major atmospheric and climatic changes."
    print(f"A. {options['A']}:\n   {analysis_a}\n")

    # Option B
    analysis_b = "Incorrect. Drought causes stomata to close, which REDUCES discrimination against 13C. This would lead to an INCREASE in the 13C ratio, which is the opposite of the observed declining trend."
    print(f"B. {options['B']}:\n   {analysis_b}\n")
    
    # Option C
    analysis_c = "Incorrect. Starch reserves are an internal energy management mechanism and reflect past conditions; they do not drive a long-term external trend in isotopic composition."
    print(f"C. {options['C']}:\n   {analysis_c}\n")

    # Option D
    analysis_d = "Incorrect. This is a secondary characteristic of tree growth, influenced by seasonal climate. It is a symptom of environmental change, not the primary cause of a centennial trend."
    print(f"D. {options['D']}:\n   {analysis_d}\n")

    # Option E
    analysis_e = "Correct. The SE Asia monsoon is a major driver of regional climate, especially rainfall. Increased monsoon intensity or duration would lead to higher water availability. This allows trees to increase stomatal conductance, leading to a decline in their 13C ratio. This large-scale climatic system integrates both regional moisture patterns and the transport of air masses carrying the global Suess effect signal, making it the predominant factor among the choices."
    print(f"E. {options['E']}:\n   {analysis_e}\n")

    # Final Conclusion
    final_answer = 'E'
    print("-" * 70)
    print(f"Final Conclusion: Based on the analysis, the most comprehensive and predominant factor is changes in the large-scale climate system.")
    print(f"\nThe final answer is: {final_answer}")

if __name__ == "__main__":
    analyze_tree_ring_data()