import textwrap

def analyze_weld_pool_flow():
    """
    Analyzes the forces in a weld pool to determine the dominant mechanism
    for inward flow in 304 stainless steel.
    """

    # --- Problem Definition ---
    material = "304 Stainless Steel"
    observation = "Outer portions of the weld pool flow inwards."
    arc_current = 350 # Amperes

    # --- Knowledge Base of Physical Forces in Welding ---
    forces = {
        'Marangoni Force': {
            'description': 'Caused by a gradient in surface tension (gamma) on the liquid surface. The flow direction depends on the temperature coefficient of surface tension (d(gamma)/dT).',
            'flow_in_pure_metal': 'Outward (d(gamma)/dT < 0)',
            'flow_in_alloy_with_surfactants': 'Inward (d(gamma)/dT > 0)'
        },
        'Arc Drag Force': {
            'description': 'Caused by the shear stress from the high-velocity plasma jet impacting the surface.',
            'flow_in_pure_metal': 'Outward (from the center)',
            'flow_in_alloy_with_surfactants': 'Outward (from the center)'
        },
        'Arc Pressure Force': {
            'description': 'Caused by the pressure of the arc plasma column pushing down on the weld pool surface, strongest at the center.',
            'flow_in_pure_metal': 'Outward (from the center)',
            'flow_in_alloy_with_surfactants': 'Outward (from the center)'
        },
        'Lorentz (electromagnetic) Force': {
            'description': 'Caused by the interaction of the welding current with its own magnetic field within the molten metal.',
            'flow_in_pure_metal': 'Inward and downward ("pinch effect")',
            'flow_in_alloy_with_surfactants': 'Inward and downward ("pinch effect")'
        },
        'Buoyancy Force': {
            'description': 'Caused by density differences due to temperature gradients. Hotter, less dense fluid rises.',
            'flow_in_pure_metal': 'Generally weak, contributing to overall circulation but not typically dominant.',
            'flow_in_alloy_with_surfactants': 'Generally weak, contributing to overall circulation but not typically dominant.'
        }
    }

    # --- Reasoning ---
    print("Step 1: Analyzing the observed phenomenon.")
    print(f"Material: {material}")
    print(f"Observation: {observation}\n")

    print("Step 2: Evaluating potential dominant forces based on the observation.\n")

    # The logic to determine the dominant force
    dominant_force = None
    reasoning = ""

    for force_name, properties in forces.items():
        flow_direction = properties['flow_in_alloy_with_surfactants']
        if "Inward" in flow_direction:
            print(f"- Considering '{force_name}': This force can cause inward flow.")
        else:
            print(f"- Considering '{force_name}': This force typically causes outward flow and does not match the observation.")

    print("\nStep 3: Concluding the dominant mechanism.\n")

    final_reasoning = (
        "Both Lorentz Force and Marangoni Force can cause inward flow. "
        "However, the composition of the alloy is critical. 304 stainless steel "
        "contains surface-active elements like sulfur. These elements reverse the "
        "typical temperature dependency of surface tension, causing surface tension to be "
        "highest in the hot center of the pool. This creates a strong surface flow from the "
        "cooler outer edges inward. This phenomenon is a well-established and powerful "
        "aspect of the Marangoni effect in alloys and is considered the dominant "
        "driver for this type of flow pattern and the resulting weld bead shape."
    )
    
    wrapper = textwrap.TextWrapper(width=80)
    print(wrapper.fill(text=final_reasoning))

    # --- Final Answer in "Equation" format ---
    # The 'numbers' are the material '304' and the property 'd(gamma)/dT > 0'
    dominant_force = 'Marangoni Force'
    factor1 = '304 Stainless Steel'
    factor2 = 'd(gamma)/dT > 0'
    
    print("\n--- Final Conclusion ---")
    print("Dominant Mechanism Equation:")
    print(f"Observed Inward Flow in '{factor1}' = {dominant_force} (driven by {factor2})")

# Execute the analysis
analyze_weld_pool_flow()