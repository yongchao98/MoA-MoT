import textwrap

def explain_sex_linked_differentiation():
    """
    Explains why high Fst between males and females points to sex chromosomes.
    """
    print("Analyzing the cause of high genetic differentiation (Fst) between males and females in a population.\n")

    explanation = {
        "A. Genetic load": "This is the reduced fitness of a population due to deleterious alleles. While it can differ between sexes, it's a consequence, not the direct cause of high differentiation at specific markers.",
        "B. XY vs ZW sex determining systems": "This is the most likely cause. Markers on sex chromosomes show extreme differentiation:\n" \
                                          "\t- In XY systems, markers on the Y-chromosome are unique to males.\n" \
                                          "\t- In ZW systems, markers on the W-chromosome are unique to females.\n" \
                                          "\t- This uniqueness leads to an Fst value of 1, the maximum possible differentiation.\n" \
                                          "\t- Markers on X and Z chromosomes also differ in frequency due to different copy numbers between sexes.",
        "C. Reproductive isolation": "This applies to different populations that cannot interbreed. Males and females within a population are, by definition, not reproductively isolated.",
        "D. Local adaptation to thermal environments": "This relates to a population adapting to its environment. It doesn't explain such a fundamental, marker-specific genetic split between sexes within that population.",
        "E. Hybrid zone dynamics": "This involves the interbreeding of two distinct populations. The question concerns differences within a single population."
    }

    print("Evaluating the potential explanations:")
    for option, desc in explanation.items():
        print(f"\n--- {option} ---")
        # textwrap.fill makes the output cleaner in the terminal
        print(textwrap.fill(desc, width=80))

    print("\n-----------------------------------------")
    print("Conclusion: The most direct and powerful explanation is that the highly differentiated markers are located on sex chromosomes.")
    print("The correct answer choice is B.")

if __name__ == "__main__":
    explain_sex_linked_differentiation()
