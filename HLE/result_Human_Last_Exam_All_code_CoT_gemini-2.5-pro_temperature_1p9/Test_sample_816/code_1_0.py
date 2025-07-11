def analyze_weld_pool_flow(material_name, dgamma_dT):
    """
    Analyzes the dominant surface flow direction in a weld pool based on the
    Marangoni effect.

    Args:
        material_name (str): The name of the material being welded.
        dgamma_dT (float): The surface tension temperature coefficient.
                           (Change in surface tension / Change in temperature)
    """
    print(f"Analyzing weld pool for: {material_name}")
    print(f"The surface tension temperature coefficient (dγ/dT) is: {dgamma_dT}")
    
    print("\nThe Marangoni force drives fluid from areas of low surface tension to high surface tension.")
    
    # In welding, the center is hotter than the edges.
    # If dgamma_dT > 0, the hot center has higher surface tension.
    # If dgamma_dT < 0, the hot center has lower surface tension.
    
    if dgamma_dT > 0:
        flow_direction = "inwards (from cooler edges to the hot center)"
        dominant_mechanism = "Marangoni Force"
        explanation = "This is typical for alloys with surfactants, like sulfur in stainless steel."
    elif dgamma_dT < 0:
        flow_direction = "outwards (from hot center to cooler edges)"
        dominant_mechanism = "Marangoni Force"
        explanation = "This is typical for pure metals."
    else:
        flow_direction = "negligible"
        dominant_mechanism = "Other forces (e.g., Lorentz, Arc Drag) would dominate surface flow."
        explanation = "No surface tension gradient exists."

    print(f"\nPredicted Surface Flow Direction: {flow_direction}")
    print(f"Explanation: {explanation}")
    print(f"The dominant mechanism for this surface flow pattern is the: {dominant_mechanism}")

# For 304 stainless steel, surfactants like sulfur reverse the coefficient, making it positive.
# We'll use a representative positive value. There isn't a single equation with numbers
# to solve, but we can model the physical principle.
material = "304 Stainless Steel"
# A positive coefficient indicates that surface tension increases with temperature.
surface_tension_temp_coefficient = 3.5e-4 # A representative positive value [N/(m*K)]

analyze_weld_pool_flow(material, surface_tension_temp_coefficient)