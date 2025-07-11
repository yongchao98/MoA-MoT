def calculate_spps_yield():
    """
    Calculates the theoretical maximum yield for a peptide synthesized
    via Solid-Phase Peptide Synthesis (SPPS).

    This demonstrates why SPPS is not suitable for long peptides like the
    100-amino-acid one in the user's question.
    """
    peptide_length = 100
    # A 99.0% yield per coupling step is optimistic; in practice, it can be lower.
    step_yield_percent = 99.0

    # The number of coupling reactions is one less than the peptide length.
    num_coupling_steps = peptide_length - 1
    step_yield_decimal = step_yield_percent / 100.0

    # The final yield is the step yield raised to the power of the number of steps.
    final_yield = step_yield_decimal ** num_coupling_steps
    final_yield_percent = final_yield * 100

    print("--- SPPS Yield Calculation ---")
    print(f"Peptide Length: {peptide_length} amino acids")
    print(f"Assumed Yield per Step: {step_yield_percent}%")
    print(f"Number of Coupling Steps: {num_coupling_steps}")
    print("\nThe final theoretical yield is calculated as:")
    # The final print statement shows the full equation as requested.
    print(f"Final Yield = ({step_yield_decimal})^{num_coupling_steps}")
    print(f"Result: {final_yield_percent:.2f}%")
    print("\nThis low theoretical yield makes purification extremely difficult and")
    print("demonstrates why alternative methods like Native Chemical Ligation are preferred.")

calculate_spps_yield()