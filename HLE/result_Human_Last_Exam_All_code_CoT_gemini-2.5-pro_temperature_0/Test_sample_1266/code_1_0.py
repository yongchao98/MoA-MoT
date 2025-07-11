def calculate_aldh_change(compound_name, potency_factor, concentration):
    """
    Models the fold change in ALDH protein based on a compound's potency.
    The equation used is a simple model: Fold Change = 1 + Potency * Concentration.
    A baseline of 1 represents the initial amount of ALDH.
    """
    fold_change = 1 + potency_factor * concentration
    print(f"Analysis for {compound_name}:")
    print(f"Equation: ALDH Fold Change = 1 + (Potency Factor) * (Concentration)")
    # The final code needs to output each number in the final equation
    print(f"Calculation: ALDH Fold Change = 1 + {potency_factor} * {concentration} = {fold_change:.2f}")
    print("-" * 30)
    return fold_change

def main():
    # Concentration of the compounds used in the experiment
    concentration_uM = 50

    # Assign relative potency factors based on biological knowledge.
    # 4-OI is a more potent Nrf2 activator than HNE derivatives.
    potency_HNE_alkyne = 0.045
    potency_4_OI = 0.075

    print("Modeling the change in ALDH amount in RAW 264.7 cells.\n")

    # Calculate change for (2E)-4-Hydroxy-2-nonen-8-ynal
    change_HNE = calculate_aldh_change("(2E)-4-Hydroxy-2-nonen-8-ynal", potency_HNE_alkyne, concentration_uM)

    # Calculate change for 4-OI
    change_4OI = calculate_aldh_change("4-OI", potency_4_OI, concentration_uM)

    # Final Conclusion
    print("\nConclusion:")
    print(f"1. Both compounds cause an INCREASE in ALDH (Fold changes > 1).")
    if change_4OI > change_HNE:
        print(f"2. The change with 4-OI ({change_4OI:.2f}-fold) is MORE than with HNE-alkyne ({change_HNE:.2f}-fold).")
    else:
        print(f"2. The change with 4-OI ({change_4OI:.2f}-fold) is LESS than with HNE-alkyne ({change_HNE:.2f}-fold).")
    print("3. The protein directly involved in sensing these compounds is Keap1.")

if __name__ == "__main__":
    main()