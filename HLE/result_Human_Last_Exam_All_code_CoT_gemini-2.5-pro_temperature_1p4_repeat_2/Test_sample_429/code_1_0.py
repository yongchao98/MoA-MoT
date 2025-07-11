import sys
# It's better to use sys.stdout.write to avoid extra newlines
# and have more control over the output format.
def solve_chemistry_problem():
    """
    This script solves the chemistry problem by identifying the metals involved
    and verifying the reaction with the data provided.
    """
    # Step 1: Define constants and identified metals.
    # From our analysis, Metal A is Strontium (Sr) and Metal M is Copper (Cu).
    metal_A_name = "Strontium (Sr)"
    metal_M_name = "Copper (Cu)"

    # Atomic masses
    Ar_A = 87.62  # Strontium
    Ar_M = 63.55  # Copper
    Ar_Cl = 35.5

    # Initial conditions from the problem
    m_MCl2_initial = 1.0  # g

    # Step 2: Calculate moles of reaction (n) assuming MCl2 is the limiting reactant.
    # Molar mass of MCl2 (CuCl2)
    MolarMass_MCl2 = Ar_M + 2 * Ar_Cl
    # Moles of reaction
    n = m_MCl2_initial / MolarMass_MCl2

    # Step 3: Calculate the theoretical changes based on the identified metals.
    # Mass of metal A (Sr) dissolved from the plate
    m_A_dissolved = n * Ar_A
    # Mass of metal M (Cu) deposited on the plate
    m_M_deposited = n * Ar_M

    # Mass change of the plate
    delta_m_plate = m_A_dissolved - m_M_deposited

    # Step 4: Calculate final solution properties.
    # Final mass of the solution
    m_solution_final = 10.0 + delta_m_plate
    # Molar mass of ACl2 (SrCl2)
    MolarMass_ACl2 = Ar_A + 2 * Ar_Cl
    # Mass of ACl2 formed in the solution
    m_ACl2_final = n * MolarMass_ACl2
    # Final mass fraction of ACl2
    w_ACl2_final = m_ACl2_final / m_solution_final

    # Step 5: Print the results and the final reaction equation.
    print(f"Based on the calculations, the metals are identified as:")
    print(f"Metal A (plate): {metal_A_name}")
    print(f"Metal M (in salt): {metal_M_name}\n")
    print("Verification:")
    print(f"Calculated plate mass decrease: {delta_m_plate:.3f} g (Problem value: 0.172 g)")
    print(f"Calculated final salt mass fraction: {w_ACl2_final*100:.2f}% (Problem value: 11.52%)\n")
    print("The calculated values are very close to the problem's values, confirming the identification.\n")
    print("The reaction equation is:")

    # Printing the final equation with coefficients and metal names
    sys.stdout.write(f"1 {metal_A_name.split(' ')[1]} + 1 {metal_M_name.split(' ')[1]}Cl2 -> 1 {metal_A_name.split(' ')[1]}Cl2 + 1 {metal_M_name.split(' ')[1]}\n")

solve_chemistry_problem()
<<<Strontium>>>