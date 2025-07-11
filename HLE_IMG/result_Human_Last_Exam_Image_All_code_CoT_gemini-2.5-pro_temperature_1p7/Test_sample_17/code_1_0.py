import math

def calculate_stress_at_notch_tip():
    """
    Explains and demonstrates the theoretical stress at the tip of a sharp wedge.
    """
    # Define symbolic variables to represent the physical quantities
    sigma_y_str = "σ_y"        # Far-field nominal stress
    kt_str = "Kt"              # Stress concentration factor
    sigma_max_str = "σ_max"    # Maximum stress at the notch tip
    rho_str = "ρ"              # Radius of curvature at the notch tip
    d_str = "d"                # Depth of the notch

    # --- Step 1: Explain the basic concept of stress concentration ---
    print("1. The Problem: Stress Concentration")
    print("The sharp wedges on the plate act as stress concentrators.")
    print("The stress at the sharp tip (point A) will be significantly higher than the applied stress (σ_y).\n")

    # --- Step 2: Introduce the governing equation ---
    print("2. The Governing Equation")
    print(f"The maximum stress at the tip ({sigma_max_str}) is found using the stress concentration factor ({kt_str}) and the nominal applied stress ({sigma_y_str}):")
    print(f"-->  {sigma_max_str} = {kt_str} * {sigma_y_str}\n")

    # --- Step 3: Explain the Stress Concentration Factor ---
    print("3. The Stress Concentration Factor (Kt)")
    print(f"For a V-notch, {kt_str} depends on the notch geometry, primarily its depth ({d_str}) and its tip radius ({rho_str}). A common formula is:")
    print(f"-->  {kt_str} = 1 + 2 * sqrt({d_str} / {rho_str})\n")

    # --- Step 4: Analyze the 'Theoretical' Case (Perfectly Sharp Tip) ---
    print("4. The Theoretical Case: A Perfectly Sharp Tip")
    print("A theoretically perfect, sharp tip has a radius of curvature equal to zero.")
    print(f"Let's substitute the value {rho_str} = 0 into our equation for {kt_str}:\n")

    # --- Step 5: Show the final equation and conclusion ---
    print("5. Final Equation and Conclusion")
    print("The full equation for the maximum stress is:")
    final_equation_str = f"{sigma_max_str} = (1 + 2 * sqrt({d_str} / {rho_str})) * {sigma_y_str}"
    print(f"-->  {final_equation_str}\n")

    print(f"When we consider a perfectly sharp tip, we analyze the limit as {rho_str} -> 0:")
    print(f"lim({rho_str}->0) [ (1 + 2 * sqrt({d_str} / {rho_str})) * {sigma_y_str} ]\n")

    print(f"Since the term 'sqrt({d_str} / {rho_str})' approaches infinity as {rho_str} approaches 0, the entire expression goes to infinity.")
    print("\n-------------------------------------------------------------")
    print("Therefore, the theoretical stress at the tip of a perfectly sharp wedge is infinite.")
    print("-------------------------------------------------------------")
    print("\nNote: In reality, materials yield or fracture before infinite stress is reached, and tips are never perfectly sharp.")


# Execute the function to print the explanation.
calculate_stress_at_notch_tip()