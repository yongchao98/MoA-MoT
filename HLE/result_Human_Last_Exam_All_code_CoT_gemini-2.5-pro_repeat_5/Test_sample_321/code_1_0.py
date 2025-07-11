import textwrap

def solve_beetle_optics():
    """
    This function analyzes the multiple-choice question about Protaetia cuprea
    and determines the most accurate structure-ecology relationship.
    """
    # Define the core facts based on biological research on scarab beetles like Protaetia cuprea.
    fact_structure = "Bouligand structures"
    fact_optical_property = "Circular polarization of light"
    fact_ecological_function = "mate attraction"

    # The question asks for the accurate structure-ecology relationship.
    # Let's break down the reasoning.
    reasoning = [
        "1. The insect in question, Protaetia cuprea, is a scarab beetle. The cuticle of many scarab beetles is known to have a specific microscopic arrangement of chitin fibers in a helical or twisted pattern. This is known as a 'Bouligand structure'.",
        "2. This Bouligand structure acts as a selective mirror for circularly polarized light. It reflects light of a specific wavelength and a specific 'handedness' (left or right) of circular polarization.",
        "3. Therefore, the key physical phenomena are 'Bouligand structures' and 'Circular polarization'. This immediately narrows down the plausible options.",
        "4. We evaluate the options with these keywords:",
        "   - Option F suggests linear polarization, which is incorrect.",
        "   - Option I suggests attracting predator attention. While possible, it's a deflection strategy.",
        "   - Option G suggests appearing unpolarized. This is a crypsis (hiding) strategy.",
        "   - Option K suggests mate attraction. This is an intraspecific communication strategy.",
        "5. Scientific studies have demonstrated that some scarab beetles can perceive circularly polarized light, which they would not be able to do if it were just for hiding from predators (who mostly can't see it). This ability to see the signal strongly supports the hypothesis that it is used as a private communication channel for recognizing potential mates.",
        "6. Therefore, the most accurate and well-supported relationship is that the Bouligand structures create circularly polarized light for mate attraction."
    ]

    correct_option_letter = "K"
    correct_option_text = "Bouligand structures - Circular polarization of  light for mate attraction"

    print("Step-by-step reasoning:")
    for step in reasoning:
        # Use textwrap to make the output clean
        print(textwrap.fill(step, width=80))
    
    print("\n" + "="*30)
    print("          Conclusion")
    print("="*30)
    print(f"The correct structure is: {fact_structure}")
    print(f"The correct optical property is: {fact_optical_property}")
    print(f"The most supported ecological function is: {fact_ecological_function}")
    print("\nMatching these facts, the final answer is:")
    print(f"Option {correct_option_letter}: {correct_option_text}")

solve_beetle_optics()