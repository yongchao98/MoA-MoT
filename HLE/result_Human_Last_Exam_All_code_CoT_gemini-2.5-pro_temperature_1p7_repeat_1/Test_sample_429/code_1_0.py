import math

def solve_chemistry_problem():
    """
    Solves the chemistry problem by testing candidate metals against the given data.
    """
    # --- Step 1: Define initial conditions from the problem statement ---
    m_solution_initial = 10.0  # g
    w_salt_initial = 0.10      # 10%
    plate_mass_decrease_given = 0.172  # g
    w_salt_final_given = 0.1152       # 11.52%
    m_salt_initial = m_solution_initial * w_salt_initial # Mass of BCl2 is 1.0 g

    print("Analyzing the problem to find the unknown divalent metal 'A'.")
    print("Assuming the reaction is A + BCl2 -> ACl2 + B\n")

    # --- Step 2: Define a function to test candidate metals ---
    def test_metals(name_A, M_A, name_B, M_B):
        """
        Calculates the expected outcomes for a given pair of metals.
        """
        M_Cl = 35.453
        M_BCl2 = M_B + 2 * M_Cl
        M_ACl2 = M_A + 2 * M_Cl

        # Moles reacted is based on the limiting reactant, BCl2
        n_moles_reacted = m_salt_initial / M_BCl2

        # Calculate theoretical mass changes
        m_A_reacted = n_moles_reacted * M_A
        m_B_deposited = n_moles_reacted * M_B
        plate_mass_decrease_calc = m_A_reacted - m_B_deposited

        # Calculate theoretical final solution properties
        m_ACl2_formed = n_moles_reacted * M_ACl2
        m_solution_final = m_solution_initial + plate_mass_decrease_calc
        w_ACl2_final_calc = m_ACl2_formed / m_solution_final
        
        print(f"--- Testing Hypothesis: A = {name_A}, B = {name_B} ---")
        print(f"Calculated Plate Mass Decrease: {plate_mass_decrease_calc:.4f} g (Given: {plate_mass_decrease_given} g)")
        print(f"Calculated Final Mass Fraction: {w_ACl2_final_calc*100:.2f}% (Given: {w_salt_final_given*100:.2f}%)")
        
        return plate_mass_decrease_calc, w_ACl2_final_calc

    # --- Step 3: Test plausible candidate pairs using rounded molar masses ---
    # As an exact algebraic solution is indeterminate, we test which common metals fit best.
    # The two most likely candidates are Ca/Mg and Sr/Cu. We use rounded molar masses
    # as is common for this type of problem.
    print("Testing candidate metal pairs with rounded molar masses...\n")

    # Candidate Pair 1: Calcium (A) and Magnesium (B)
    ca_dm, ca_w = test_metals("Calcium", 40, "Magnesium", 24)

    # Candidate Pair 2: Strontium (A) and Copper (B)
    sr_dm, sr_w = test_metals("Strontium", 88, "Copper", 64)
    
    # --- Step 4: Compare results and conclude ---
    # Calculate the error for each pair to see which is a better fit.
    error_dm_ca = abs(ca_dm - plate_mass_decrease_given)
    error_w_ca = abs(ca_w - w_salt_final_given)

    error_dm_sr = abs(sr_dm - plate_mass_decrease_given)
    error_w_sr = abs(sr_w - w_salt_final_given)

    print("\n--- Conclusion ---")
    if (error_dm_ca + error_w_ca) < (error_dm_sr + error_w_sr):
        metal_A_name = "Calcium"
        metal_A_symbol = "Ca"
        metal_B_symbol = "Mg"
        print("The data best fits the reaction between Calcium and Magnesium Chloride.")
    else:
        metal_A_name = "Strontium"
        metal_A_symbol = "Sr"
        metal_B_symbol = "Cu"
        print("The data best fits the reaction between Strontium and Copper Chloride.")
        
    print(f"\nThe unknown metal (A) is: {metal_A_name}")
    
    # --- Step 5: Print the final chemical equation with coefficients ---
    print("\nThe equation for the reaction is:")
    print(f"1{metal_A_symbol} + 1{metal_B_symbol}Cl2 -> 1{metal_A_symbol}Cl2 + 1{metal_B_symbol}")

# Execute the function to solve the problem
solve_chemistry_problem()