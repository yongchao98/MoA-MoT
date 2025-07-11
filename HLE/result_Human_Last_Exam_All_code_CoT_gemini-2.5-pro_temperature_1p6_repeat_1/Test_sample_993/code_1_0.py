import textwrap

def explain_paleoclimate_isotope():
    """
    Explains the reasoning behind the relationship between climate and 18O isotope levels in lake sulfate deposits.
    """

    print("Step 1: Understand the Goal")
    print("The goal is to find the climate conditions that result in low levels of the heavy oxygen isotope, 18O, in lake sulfate deposits.")
    print("-" * 20)

    print("Step 2: How Oxygen Isotopes Work in Water")
    explanation_isotopes = """
    The oxygen in lake sulfate (SO4) comes from the lake's water (H2O). Therefore, the 18O level in the sulfate reflects the 18O level in the water.
    - Evaporation: Water with the lighter 16O isotope evaporates more easily. This leaves the remaining lake water enriched in the heavier 18O isotope. High evaporation occurs in WARM and DRY climates.
    - Precipitation: Water vapor loses heavy 18O as it cools and travels. Therefore, precipitation, especially in COLD climates, is depleted in 18O (has low 18O levels).
    """
    print(textwrap.dedent(explanation_isotopes).strip())
    print("-" * 20)

    print("Step 3: Combine the Factors for LOW 18O Levels")
    explanation_synthesis = """
    To get LOW 18O levels in the sulfate, we need LOW 18O levels in the lake water. This requires:
    - Minimizing evaporation: This points to a COLD climate.
    - Maximizing input of 18O-depleted water: This points to a WET climate (high precipitation).
    - Lake Levels: HIGH lake levels are a consequence of a wet climate where input exceeds output (evaporation). A deep lake is also less affected by evaporation than a shallow one.
    """
    print(textwrap.dedent(explanation_synthesis).strip())
    print("-" * 20)

    print("Step 4: Conclusion")
    print("The combination of a Wet, Cold Climate with High Lake Levels creates the ideal environment for producing sulfate deposits with low levels of the 18O isotope.")
    print("This matches option G.")
    print("-" * 20)

    print("Final Answer:")
    print("G")

# Run the explanation
explain_paleoclimate_isotope()