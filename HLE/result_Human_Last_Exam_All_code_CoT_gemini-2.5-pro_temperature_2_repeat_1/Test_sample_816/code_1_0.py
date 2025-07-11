import textwrap

def analyze_welding_flow():
    """
    Analyzes the dominant mechanism for inward weld pool flow in GTAW
    of 304 stainless steel.
    """
    # Problem statement
    observation = "The outer portions of the weld pool flow inwards from the cooler edge to the hotter center."
    material = "304 stainless steel"
    process = "Gas Tungsten Arc Welding (GTAW)"

    # Answer choices
    choices = {
        'A': "Marangoni Force",
        'B': "Arc drag Force",
        'C': "Arc Pressure Force",
        'D': "Lorentz (electromagnetic) Force",
        'E': "Buoyancy Force"
    }

    # Explanation
    print(f"Process: {process}")
    print(f"Material: {material}")
    print(f"Observation: {observation}\n")
    print("--- Analysis of Driving Forces ---")

    explanation = f"""
    1. The direction of fluid flow on the weld pool surface is determined by the dominant force acting upon it. The key is to find the force that causes flow from the cooler edge TOWARDS the hot center.

    2. The Marangoni Force arises from a gradient in surface tension (γ). Fluid flows from a region of low surface tension to a region of high surface tension. The relationship between surface tension and temperature (T) is critical.

       - For most pure metals, surface tension decreases with temperature (dγ/dT < 0). This would cause an OUTWARD flow from the hot center to the cool edge.

       - However, {material} contains surface-active elements like sulfur. These elements cause the surface tension to INCREASE with temperature (dγ/dT > 0) in the typical welding temperature range.

       - A positive dγ/dT means surface tension is highest at the hot center and lowest at the cooler edge. This drives the observed INWARD flow.

    3. Other forces do not fit the observation as the primary cause:
       - {choices['B']}: Pushes fluid OUTWARD.
       - {choices['D']}: Creates an inward and downward 'pinch' but the Marangoni effect dictates the surface flow direction.
       - {choices['C']} & {choices['E']}: Are considered less dominant than the Marangoni force for determining the surface flow pattern in this scenario.

    Conclusion: The observed inward flow is a classic signature of the Marangoni effect in materials with a positive temperature coefficient of surface tension.
    """

    # Print formatted explanation and the final answer
    print(textwrap.dedent(explanation))
    
    final_answer_key = 'A'
    final_answer_text = choices[final_answer_key]
    print(f"Therefore, the dominant mechanism is: {final_answer_key}. {final_answer_text}")

# Run the analysis
analyze_welding_flow()