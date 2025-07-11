import textwrap

def explain_welding_flow():
    """
    Analyzes the forces in a weld pool to determine the dominant mechanism
    for the observed inward flow in stainless steel.
    """

    # Define the core problem from the user's query
    observation = "The outer portions of the weld pool flow inwards."
    material = "304 stainless steel"
    process = "Gas Tungsten Arc Welding (GTAW)"

    # Print a structured explanation
    print(f"Analysis of weld pool fluid dynamics for {process} of {material}:")
    print("-" * 70)
    print(f"Key Observation: {observation}")
    print("-" * 70)

    # Explanation of forces and their typical effect
    explanation = {
        'A. Marangoni Force': 'Flow driven by surface tension gradients. The direction depends on material composition. In steels with surfactants (like sulfur), it causes INWARD flow from the cooler edges to the hot center.',
        'B. Arc drag Force': 'Shear force from the plasma gas jet. It pushes the liquid OUTWARD from the center.',
        'C. Arc Pressure Force': 'Pressure from the plasma jet. It depresses the pool center and pushes liquid OUTWARD.',
        'D. Lorentz (electromagnetic) Force': 'Interaction between current and magnetic field (pinch effect). It drives the bulk liquid INWARD and downward.',
        'E. Buoyancy Force': 'Flow driven by density differences due to heat. Generally a weaker effect compared to Marangoni or Lorentz forces.'
    }

    print("Evaluating the potential forces:\n")
    for force, desc in explanation.items():
        # Use textwrap for clean formatting
        wrapped_desc = textwrap.fill(desc, width=68, initial_indent='    ', subsequent_indent='    ')
        print(f"  - {force}:\n{wrapped_desc}\n")

    print("-" * 70)
    print("Conclusion:")
    print("-" * 70)

    conclusion_text = (
        "The observed INWARD flow from the outer portions is a surface phenomenon. "
        "While the Lorentz force causes bulk inward flow, the Marangoni force is "
        "the primary driver of this specific surface flow pattern in stainless steels.\n\n"
        "This happens because trace elements like sulfur, present in 304 stainless steel, "
        "act as surfactants. They reverse the normal behavior of surface tension, causing it "
        "to be highest at the hottest point (the center). This high surface tension in the "
        "center pulls the surface liquid from the cooler edges inward.\n\n"
        "Therefore, the dominant mechanism for the observed flow is the Marangoni Force."
    )
    print(textwrap.fill(conclusion_text, width=70))
    print("\n")


# Execute the explanation function
explain_welding_flow()
<<<A>>>