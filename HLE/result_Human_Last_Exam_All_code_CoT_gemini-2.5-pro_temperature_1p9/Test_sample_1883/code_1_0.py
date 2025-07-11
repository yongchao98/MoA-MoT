import math

def predict_product_ratio():
    """
    Predicts the product ratio of a thermal electrocyclic reaction
    based on Frontier Molecular Orbital theory and conformer analysis.
    """

    # Step 1 & 2: Theoretical Background (as explained in the plan)
    # The reaction is an 8-pi electron thermal electrocyclization, which proceeds via
    # a conrotatory mechanism according to FMO theory.
    # The ratio of cis (A) and trans (B) products depends on the relative stabilities
    # of the four possible reactive conformers of the starting material.
    
    # Step 3: Analysis of Conformers and Products
    # Let 'in' and 'out' describe the orientation of the terminal methyl groups.
    # Conrotatory motion applied to the four conformers gives:
    # 1. C2(Me-in)/C9(Me-out) -> trans-isomer (B)
    # 2. C2(Me-out)/C9(Me-out) -> cis-isomer (A)
    # 3. C2(Me-in)/C9(Me-in) -> cis-isomer (A)
    # 4. C2(Me-out)/C9(Me-in) -> trans-isomer (B)
    
    # Step 4: Estimating Conformer Stabilities
    # Stability depends on steric clashes (A(1,3) strain). An 'in' methyl group
    # clashes with an adjacent backbone hydrogen. An 'out' methyl does not.
    # E_clash: Energy cost for one 'in' methyl group. We use a standard value.
    E_clash_cal_per_mol = 1000.0  # cal/mol (equivalent to 1.0 kcal/mol)
    
    # Relative energies of the four conformers (ΔE):
    # Conf 2 (out/out): 0 clashes (most stable). ΔE = 0
    # Conf 1 (in/out) & 4 (out/in): 1 clash. ΔE = E_clash
    # Conf 3 (in/in): 2 clashes. ΔE = 2 * E_clash
    
    # Physical constants
    R_cal_per_mol_K = 1.987  # Gas constant in cal/(mol*K)
    T_K = 298.0              # Room temperature in Kelvin
    
    # Step 5: Calculating the product ratio
    # Ratio [A]/[B] = (Pop(Conf 2) + Pop(Conf 3)) / (Pop(Conf 1) + Pop(Conf 4))
    # Pop(i) is proportional to exp(-ΔE_i / RT)
    # Ratio [A]/[B] = (exp(0) + exp(-2*E_clash/RT)) / (exp(-E_clash/RT) + exp(-E_clash/RT))
    # This simplifies to cosh(E_clash / RT)
    
    # Calculate x = E_clash / (R * T)
    x = E_clash_cal_per_mol / (R_cal_per_mol_K * T_K)
    
    # Calculate the ratio A/B
    ratio_A_to_B = math.cosh(x)

    # Print the explanation and the result, showing the numbers in the equation
    print("The electrocyclic reaction of (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene is conrotatory under thermal conditions.")
    print("The ratio of the cis-isomer (A) to the trans-isomer (B) is determined by the relative stabilities of the reacting conformers.")
    print("This ratio can be calculated using the formula:")
    print("Ratio [A]/[B] = cosh(E_clash / (R * T))")
    print("\nUsing the following values:")
    print(f"Steric clash energy (E_clash): {E_clash_cal_per_mol} cal/mol")
    print(f"Gas constant (R): {R_cal_per_mol_K} cal/mol*K")
    print(f"Temperature (T): {T_K} K")
    print("\nThe calculation is:")
    print(f"Ratio [A]/[B] = cosh({E_clash_cal_per_mol} / ({R_cal_per_mol_K} * {T_K}))")
    print(f"Ratio [A]/[B] = cosh({x:.4f})")
    print(f"Ratio [A]/[B] = {ratio_A_to_B:.2f}")
    print("\nThis means the ratio of isomer A to isomer B is approximately 2.80 : 1.")
    
predict_product_ratio()