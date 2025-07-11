def explain_marangoni_flow(material_name, temp_coeff_surface_tension_sign):
    """
    Explains the direction of Marangoni flow in a weld pool
    based on the sign of the temperature coefficient of surface tension (d(gamma)/dT).
    """
    print(f"Analyzing weld pool surface flow for: {material_name}")
    print(f"The temperature coefficient of surface tension (d(gamma)/dT) is {temp_coeff_surface_tension_sign}.")

    # The center of the weld pool is hotter than the edge.
    # Flow is directed from regions of low surface tension to high surface tension.
    if temp_coeff_surface_tension_sign == "Negative":
        # d(gamma)/dT < 0 means surface tension is highest at the cool edge.
        flow_direction = "Outward (from hot center to cool edge)"
    elif temp_coeff_surface_tension_sign == "Positive":
        # d(gamma)/dT > 0 means surface tension is highest at the hot center.
        flow_direction = "Inward (from cool edge to hot center)"
    else:
        flow_direction = "Undefined"

    print(f"Resulting dominant surface flow direction: {flow_direction}")
    print("-" * 20)


# Case 1: Pure metal (hypothetical)
explain_marangoni_flow("Pure Metal", "Negative")

# Case 2: 304 Stainless Steel (as described in the problem)
# The presence of surfactants like sulfur and oxygen makes the coefficient positive.
explain_marangoni_flow("304 Stainless Steel", "Positive")

print("Conclusion: The observed inward flow in 304 stainless steel is characteristic of a positive")
print("temperature coefficient of surface tension. This phenomenon is driven by the Marangoni Force.")
