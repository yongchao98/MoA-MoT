def calculate_theory_anomaly(theory_description):
    """
    Calculates the 't Hooft anomaly coefficient for a given theory.
    This is a simplified model for illustrative purposes.
    The anomaly for the SU(N_f)_L x SU(N_f)_R global symmetry in QCD is proportional to N_c.
    """
    if theory_description.get("level") == "UV":
        # In the UV theory (quarks & gluons), the anomaly is proportional to the number of colors.
        n_c = theory_description.get("colors", 0)
        print(f"UV Theory (Fundamental): Based on {theory_description.get('quarks')} quarks and {n_c} colors.")
        # The anomaly coefficient 'A' is proportional to N_c. For simplicity, we set it to N_c.
        return n_c
    elif theory_description.get("level") == "IR":
        # In the IR theory (composite hadrons), the anomaly is matched by the
        # Wess-Zumino-Witten (WZW) term for the Goldstone bosons (pions).
        # The coefficient of this term must match the UV anomaly.
        wzw_coefficient = theory_description.get("wzw_coefficient", 0)
        print(f"IR Theory (Effective): Based on {theory_description.get('pions')} composite pions.")
        print("The anomaly is reproduced by a Wess-Zumino-Witten term.")
        return wzw_coefficient
    else:
        return 0

def main():
    """
    Main function to demonstrate the 't Hooft anomaly matching condition.
    """
    # Define the parameters for our example theory (like QCD)
    num_colors = 3
    num_flavors = 2

    # 1. Define the high-energy (UV) theory
    uv_theory = {
        "level": "UV",
        "quarks": num_flavors,
        "colors": num_colors
    }

    # 2. Define the low-energy (IR) effective theory
    # The number of pions (Goldstone bosons) is N_f^2 - 1
    num_pions = num_flavors**2 - 1
    # For the IR theory to be valid, its anomaly must match the UV theory's anomaly.
    # The WZW term coefficient must be equal to N_c.
    ir_theory = {
        "level": "IR",
        "pions": num_pions,
        "wzw_coefficient": num_colors
    }

    # 3. Calculate the anomaly in both theories
    anomaly_uv = calculate_theory_anomaly(uv_theory)
    print("-" * 30)
    anomaly_ir = calculate_theory_anomaly(ir_theory)
    print("-" * 30)

    # 4. State the matching condition and its implication
    print("The 't Hooft Anomaly Matching Condition states: Anomaly(UV) = Anomaly(IR)")
    print("This equality is a powerful non-perturbative test.")
    print("\nPhysical Implication:")
    print("Any proposed low-energy effective theory MUST reproduce the anomaly of the fundamental UV theory.")
    print("This provides a strong CONSTRAINT on the possible low-energy degrees of freedom and their interactions.")
    print("\nChecking the constraint for our example:")
    
    # Final equation with numbers
    print(f"Final Equation: {anomaly_uv} = {anomaly_ir}")

    if anomaly_uv == anomaly_ir:
        print("The condition is satisfied. The proposed IR theory is consistent.")
    else:
        print("The condition is NOT satisfied. The proposed IR theory is invalid.")

if __name__ == "__main__":
    main()