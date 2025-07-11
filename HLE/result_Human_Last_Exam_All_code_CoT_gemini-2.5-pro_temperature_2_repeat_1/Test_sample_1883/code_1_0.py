import math

def predict_product_ratio():
    """
    Calculates the ratio of cis (A) to trans (B) isomers in the thermal
    electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene.

    The prediction is based on a competition between an orbitally "allowed"
    pathway leading to the sterically hindered cis-isomer (A) and an
    orbitally "forbidden" pathway leading to the sterically favored
    trans-isomer (B).
    """

    # --- Step 1: Define physical constants and energy estimates ---

    # Estimated orbital penalty for a forbidden 8-pi electron electrocyclization (kcal/mol)
    delta_E_orbital = 15.0
    
    # Estimated steric destabilization of the "allowed" transition state leading
    # to the cis-product relative to the "forbidden" TS leading to the trans-product (kcal/mol).
    # This value is chosen to be high due to the reactant's constrained helical shape.
    delta_E_steric = 18.0
    
    # Gas constant in kcal/mol.K
    R_kcal_per_mol_K = 0.001987
    
    # Assumed reaction temperature in Kelvin (100 degrees Celsius)
    T_kelvin = 373.15

    # --- Step 2: Calculate the difference in activation energies (ΔΔG‡) ---
    
    # Ea(A) is the activation energy for the allowed path to cis-product A.
    # Ea(B) is the activation energy for the forbidden path to trans-product B.
    # We can model the difference as:
    # ΔΔG‡ = Ea(A) - Ea(B)
    #      = (Ea_allowed + ΔE_steric) - (Ea_forbidden)
    #      = ΔE_steric - ΔE_orbital
    
    delta_delta_G_kcal = delta_E_steric - delta_E_orbital
    
    # --- Step 3: Calculate the product ratio using the Boltzmann distribution ---
    # Ratio A/B = k_A / k_B = exp( - (Ea(A) - Ea(B)) / RT )
    
    exponent = -delta_delta_G_kcal / (R_kcal_per_mol_K * T_kelvin)
    ratio_A_to_B = math.exp(exponent)

    # --- Step 4: Print the detailed calculation and result ---
    
    print("Prediction of the product ratio A (cis) / B (trans):")
    print("-" * 50)
    print("The ratio is calculated from the difference in activation energies (ΔΔG‡) for the two competing pathways:")
    print("Ratio A/B = exp( - (Ea(A) - Ea(B)) / (R * T) )")
    print("\nComponent values used in the calculation:")
    print(f"  - Orbital penalty (Ea_forbidden - Ea_allowed): {delta_E_orbital} kcal/mol")
    print(f"  - Steric penalty (E_steric_A - E_steric_B): {delta_E_steric} kcal/mol")
    print(f"  - Gas constant (R): {R_kcal_per_mol_K} kcal/mol.K")
    print(f"  - Temperature (T): {T_kelvin} K")
    
    print("\nStep-by-step calculation:")
    print(f"  1. Difference in activation energies, ΔΔG‡ = Ea(A) - Ea(B)")
    print(f"     ΔΔG‡ = ΔE_steric - ΔE_orbital = {delta_E_steric} - {delta_E_orbital} = {delta_delta_G_kcal:.1f} kcal/mol")
    
    print(f"\n  2. Final Ratio Calculation:")
    final_equation = f"Ratio A/B = exp( -({delta_delta_G_kcal:.1f}) / ({R_kcal_per_mol_K} * {T_kelvin}) )"
    print(f"     {final_equation}")
    
    print(f"     Ratio A/B = exp({exponent:.4f})")
    
    print("\nFinal Result:")
    print(f"The predicted ratio of A (cis) to B (trans) is approximately 1 : {1/ratio_A_to_B:.1f}")
    
    # Return the numerical ratio for the final answer block
    return ratio_A_to_B

# Run the prediction
final_ratio = predict_product_ratio()
# The final output requires a single number for the A/B ratio.
# final_ratio represents A/B, so it is the answer.
# Let's print the value as requested.
# <<<RATIO_A/B>>>
print(f"\nValue of Ratio A/B = {final_ratio:.4f}")

if __name__ == '__main__':
    # The script can be run directly if needed
    pass

# For the final answer block, we output the numerical value of the A/B ratio
# print(f'<<<{final_ratio:.4f}>>>') is not standard format, the triple chevron needs to be the final thing.
# The calculation shows A is the minor product and B is major. 
# ratio A/B = 0.0175. So B:A is ~57:1. The ratio A to B is 0.0175
# The user wants "the ratio", which can be interpreted in a few ways.
# Let's provide the number A/B.
final_answer_value = 1/57.1
# Let's use the actual calculated value
final_answer_value = 0.01746
