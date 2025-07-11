import math

def calculate_marangoni_effect():
    """
    This function illustrates the principle of the Marangoni force in a weld pool.
    The direction of the Marangoni force is determined by the sign of the temperature
    coefficient of surface tension (dγ/dT).
    """

    print("Analyzing the dominant force for inward weld pool flow.")
    print("-" * 50)
    print("The primary force governing surface flow in a weld pool is the Marangoni force.")
    print("It is described by the Marangoni shear stress (τ) equation:")
    print("\nτ = - (dγ/dT) * ∇T\n")
    print("Where:")
    print("  τ   = Marangoni shear stress (force per unit area) along the surface.")
    print("  dγ/dT = Temperature coefficient of surface tension.")
    print("  ∇T  = Temperature gradient from the center to the edge.\n")

    # For 304 stainless steel, which contains sulfur, dγ/dT can be positive.
    # This is a known phenomenon where surface-active elements reverse the usual effect.
    d_gamma_d_T = 0.0003  # N/(m·K) (A plausible positive value for steel with surfactants)

    # Temperature is highest at the center and decreases towards the edge.
    # So, the temperature gradient (∇T = T_edge - T_center / distance) is negative.
    temperature_gradient = -100000  # K/m (An illustrative value for a steep gradient)

    print("Step 1: Define variables for 304 Stainless Steel")
    print(f"   - Let's assume a positive dγ/dT for this material: {d_gamma_d_T} N/(m·K)")
    print(f"   - The temperature gradient (∇T) from the hot center to the cool edge is negative: {temperature_gradient} K/m\n")

    # Calculate the Marangoni stress
    marangoni_stress = -1 * (d_gamma_d_T) * (temperature_gradient)

    print("Step 2: Calculate the Marangoni stress using the final equation")
    print(f"τ = - ({d_gamma_d_T}) * ({temperature_gradient})")
    print(f"τ = {marangoni_stress:.2f} N/m^2 (Pascals)\n")

    print("Step 3: Interpret the result")
    if marangoni_stress > 0:
        print("The positive result for the stress (τ) indicates the force is directed against the temperature gradient.")
        print("This means the flow is from the cooler edge towards the hotter center.")
        print("This 'inward' flow is exactly what is described in the problem.")
    else:
        print("A negative result would indicate a force in the direction of the gradient (outward flow).")
    
    print("\nConclusion: The Marangoni force is the dominant mechanism for this observed inward flow.")

calculate_marangoni_effect()