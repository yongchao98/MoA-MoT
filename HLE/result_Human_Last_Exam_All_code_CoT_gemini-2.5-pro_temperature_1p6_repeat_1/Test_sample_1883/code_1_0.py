import math

def predict_isomer_ratio():
    """
    Predicts the ratio of cis (A) to trans (B) isomers in the thermal
    electrocyclization of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene using
    Frontier Molecular Orbital theory and steric considerations.
    """

    # --- Define Constants ---

    # The activation energy difference (dE) between the transition states
    # leading to the cis and trans products. The cis transition state is higher
    # in energy due to a methyl-methyl steric clash, estimated at 0.9 kcal/mol.
    # We convert it to cal/mol for consistency with the gas constant R.
    delta_E_kcal_per_mol = 0.9
    delta_E_cal_per_mol = delta_E_kcal_per_mol * 1000

    # The ideal gas constant (R) in cal/mol-K.
    R_cal_per_mol_K = 1.987

    # The absolute temperature (T) in Kelvin. A literature temperature
    # for this reaction under "thermal conditions" is 110 degrees Celsius.
    T_celsius = 110.0
    T_kelvin = T_celsius + 273.15

    # --- Calculation ---

    # The ratio of product A to product B is given by the equation:
    # Ratio = exp(-dE / (R * T))
    exponent_value = -delta_E_cal_per_mol / (R_cal_per_mol_K * T_kelvin)
    ratio_A_to_B = math.exp(exponent_value)

    # --- Output Results ---

    print("Prediction of the product ratio (A/B) based on FMO theory and transition state energies:")
    print("-" * 75)
    print("The final equation for the ratio is: A/B = exp(-dE / (R * T))")
    print("\nUsing the following values:")
    print(f"  dE (Energy Difference) = {delta_E_cal_per_mol} cal/mol")
    print(f"  R (Gas Constant)       = {R_cal_per_mol_K} cal/mol-K")
    print(f"  T (Temperature)        = {T_kelvin:.2f} K ({T_celsius}°C)")
    
    print("\nSubstituting the values into the equation:")
    print(f"  A/B = exp(-{delta_E_cal_per_mol} / ({R_cal_per_mol_K} * {T_kelvin:.2f}))")
    print(f"  A/B = exp({exponent_value:.4f})")
    
    print("\n--- Final Predicted Ratio ---")
    print(f"The predicted ratio of cis-isomer A to trans-isomer B is approximately: {ratio_A_to_B:.3f}")
    
    if ratio_A_to_B > 0:
        ratio_B_to_A = 1 / ratio_A_to_B
        print(f"This corresponds to a B:A ratio of approximately {ratio_B_to_A:.1f} : 1.")

if __name__ == '__main__':
    predict_isomer_ratio()