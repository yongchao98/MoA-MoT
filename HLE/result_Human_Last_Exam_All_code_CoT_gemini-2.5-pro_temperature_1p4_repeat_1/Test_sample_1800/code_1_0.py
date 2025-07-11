import textwrap

def get_catalyst_info():
    """
    Prints information on the ideal Ni/Ce ratio for catalytic reactions
    based on a review of scientific literature.
    """
    
    print("Analysis of Ideal Ni/Ce Ratio for Ni-Ceria Nanoparticle Catalysts")
    print("="*70)
    
    # General Context
    intro_text = """
    The ideal Ni/Ce ratio in Ni-Ceria (Ni-CeO2) catalysts is not a single, universally fixed value. It depends on synthesis methods, reaction conditions (temperature, pressure), and particle morphology. However, a general optimal range has been identified to maximize the crucial Ni-Ceria interface, which is believed to host the active sites for the reactions, while preventing the sintering of Ni particles.
    """
    print(textwrap.dedent(intro_text))
    
    # Specific Reactions
    print("--- Ideal Ratio for Water Gas Shift (WGS) and Water Splitting (WS) ---")
    
    wgs_ws_text = """
    For both the Water Gas Shift (WGS) reaction and water splitting, a relatively low nickel loading is preferred. This ensures high dispersion of Ni atoms on the ceria support and maximizes the perimeter contact between Ni and Ceria.
    
    Research consistently shows high catalytic activity when the atomic ratio of Ni to Ce is approximately 1 to 9.
    """
    print(textwrap.dedent(wgs_ws_text))

    # The final equation as requested
    print("The ideal atomic ratio can be expressed as:")
    ni_val = 1
    ce_val = 9
    print(f"    Ni : Ce = {ni_val} : {ce_val}")
    
    # Expressing as atomic percentage
    total_parts = ni_val + ce_val
    ni_atomic_percent = (ni_val / total_parts) * 100
    print(f"\nThis corresponds to a Nickel atomic loading of Ni/(Ni+Ce) = {ni_val}/{total_parts}, or {ni_atomic_percent:.1f}%.")
    
    conclusion_text = """
    This loading balances the availability of active Ni sites with the prevention of bulk NiO formation, which is less active and prone to deactivation. While ratios up to 1:4 (20% Ni) can also be effective, the 1:9 ratio is a very common and highly cited optimum.
    """
    print(textwrap.dedent(conclusion_text))

# Run the function to display the information
if __name__ == "__main__":
    get_catalyst_info()