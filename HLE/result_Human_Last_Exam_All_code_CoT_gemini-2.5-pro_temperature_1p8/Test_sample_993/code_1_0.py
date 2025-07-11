import textwrap

def explain_paleoclimate_isotope_question():
    """
    Analyzes the factors affecting 18O isotope levels in ancient lake sulfates
    and prints a step-by-step explanation.
    """
    
    # Define the problem and the goal
    title = "Analysis of Factors for Low ¹⁸O Levels in Lake Sulfates"
    print(title)
    print("=" * len(title))
    print(textwrap.fill(
        "The goal is to identify the climate and lake conditions that result in "
        "sulfate deposits with low concentrations of the heavy oxygen isotope, ¹⁸O.",
        width=80
    ))
    print("\n--- Key Principle: Isotopic Fractionation via Evaporation ---")
    print(textwrap.fill(
        "During evaporation, water with the lighter ¹⁶O isotope evaporates more "
        "easily than water with the heavier ¹⁸O isotope. This leaves the remaining "
        "lake water enriched in ¹⁸O. Therefore, high evaporation leads to high ¹⁸O levels "
        "in the lake and its mineral deposits.",
        width=80
    ))
    
    print("\n--- Evaluating the Environmental Factors for LOW ¹⁸O ---")
    
    # Step 1: Wet vs. Dry Climate
    print("\n1. Climate: Wet vs. Dry")
    print(textwrap.fill(
        "A 'WET' climate (high precipitation/inflow) replenishes the lake with "
        "isotopically light rainwater, counteracting the effects of evaporation and keeping "
        "¹⁸O levels low. A 'DRY' climate does the opposite.",
        initial_indent="   ", subsequent_indent="   ", width=80
    ))
    
    # Step 2: Cold vs. Warm Climate
    print("\n2. Temperature: Cold vs. Warm")
    print(textwrap.fill(
        "A 'COLD' climate reduces the rate of evaporation. Less evaporation means "
        "less enrichment of ¹⁸O in the lake water. A 'WARM' climate accelerates "
        "evaporation.",
        initial_indent="   ", subsequent_indent="   ", width=80
    ))
    
    # Step 3: High vs. Shallow Lake Levels
    print("\n3. Lake Level: High vs. Shallow")
    print(textwrap.fill(
        "'HIGH' lake levels indicate a large volume of water, which is more resistant "
        "to isotopic changes from evaporation. High levels also suggest that inflow exceeds "
        "evaporation, which is consistent with a wet climate. 'SHALLOW' lakes have a "
        "high surface-area-to-volume ratio, promoting rapid evaporation.",
        initial_indent="   ", subsequent_indent="   ", width=80
    ))
    
    # Conclusion
    print("\n--- Conclusion ---")
    print(textwrap.fill(
        "To achieve LOW levels of the ¹⁸O isotope, we need conditions that minimize "
        "evaporation and maximize inflow. The ideal combination is a WET, COLD climate "
        "with HIGH lake levels.",
        width=80
    ))
    
    # Match to the answer choice
    print("\nMatching this conclusion with the given options leads to:")
    print("G. Wet, cold climate with high lake levels")

# Run the explanation function
explain_paleoclimate_isotope_question()