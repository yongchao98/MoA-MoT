import textwrap

def analyze_tree_ring_factors():
    """
    Analyzes the potential factors for a declining 13C ratio in Chinese pine tree rings
    from 1886-1990 and identifies the most predominant factor among the given choices.
    """

    # The central observation to be explained
    observation = "A declining 13C isotope ratio in tree rings over the period 1886-1990 AD."

    # The provided answer choices
    options = {
        'A': "An increase in tree ring thickness as the tree matures",
        'B': "Periods of drought",
        'C': "Increased photosynthetic reserves of starch fueling tree growth",
        'D': "Thinning earlywood tree ring proportion",
        'E': "Changes in the SE Asia monsoon"
    }

    print(f"Analyzing the Predominant Factor for: {observation}\n")
    print("-" * 70)

    # Dictionary to hold the reasoning for each option
    analysis = {
        'A': "This is a tree-level physiological effect ('ontogenetic trend') and is unlikely to be the predominant driver of a consistent, regional, century-long trend.",
        'B': "Drought causes stomatal closure to conserve water. This leads to *less* discrimination against the heavy 13C isotope, causing the 13C ratio to INCREASE, which is the opposite of the observed trend.",
        'C': "The use of starch reserves is a short-term process for managing energy and does not explain a systematic, multi-decadal decline in the isotope ratio.",
        'D': "Changes in earlywood/latewood proportions are a response to seasonal conditions, not a primary driver of a long-term trend. This is a secondary effect.",
        'E': "The SE Asia monsoon is a major regional climate system. A long-term strengthening of the monsoon leads to wetter and more humid conditions. This allows trees to keep their stomata wide open, maximizing discrimination against the heavy 13C isotope and causing a DECLINE in the 13C ratio in their wood. This matches the observation."
    }

    # Print the analysis for each option
    for choice, text in options.items():
        print(f"Choice {choice}: {text}")
        # Use textwrap for clean printing of longer reasoning strings
        wrapped_analysis = textwrap.fill(f"Analysis: {analysis[choice]}", width=70)
        print(wrapped_analysis)
        print("-" * 70)

    # Conclude with the most plausible answer
    final_answer = 'E'
    print("\nConclusion:")
    print("The most plausible predominant factor is E. A stronger monsoon creates favorable (wetter) conditions,")
    print("leading to a physiological response in the tree that results in a declining 13C ratio.")
    print(f"\nThe final answer is option: {final_answer}")

if __name__ == "__main__":
    analyze_tree_ring_factors()