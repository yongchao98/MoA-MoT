def model_bond_selective_chemistry():
    """
    This script models how exciting a C-H bond with a laser affects its
    reactivity compared to C-D bonds in a reaction with Fluorine.
    """

    # 1. Define base reactivity values.
    # Due to the kinetic isotope effect, C-H bonds are naturally
    # slightly more reactive than C-D bonds. We'll use arbitrary units.
    base_reactivity_CH = 1.5
    base_reactivity_CD = 1.0

    # 2. Define the effect of laser excitation.
    # Vibrational excitation provides a lot of energy, making the bond
    # much more likely to react. This factor is illustrative.
    excitation_enhancement_factor = 50.0

    print("--- Scenario 1: No Laser Excitation (Thermal Reaction) ---")
    # In CHD3, there is 1 C-H bond and 3 C-D bonds.
    prob_H_removal = 1 * base_reactivity_CH
    prob_D_removal = 3 * base_reactivity_CD
    total_reactivity = prob_H_removal + prob_D_removal
    print(f"Relative chance of H removal: {prob_H_removal}")
    print(f"Relative chance of D removal: {prob_D_removal}")
    print(f"Likelihood of H removal is {prob_H_removal / prob_D_removal:.2f} times that of D removal.\n")


    print("--- Scenario 2: C-H Bond Excited by Laser ---")
    # The laser energy is targeted only at the C-H bond.
    excited_reactivity_CH = base_reactivity_CH * excitation_enhancement_factor

    # We can write the calculation for the new C-H reactivity as an equation:
    print("Equation for excited C-H reactivity: New = Base * Factor")
    print(f"Output of equation: {excited_reactivity_CH:.1f} = {base_reactivity_CH:.1f} * {excitation_enhancement_factor:.1f}")

    # The reactivity of the C-D bonds is unchanged.
    prob_H_removal_excited = 1 * excited_reactivity_CH
    prob_D_removal_excited = 3 * base_reactivity_CD
    total_reactivity_excited = prob_H_removal_excited + prob_D_removal_excited

    print(f"\nRelative chance of H removal (excited): {prob_H_removal_excited:.1f}")
    print(f"Relative chance of D removal (not excited): {prob_D_removal_excited:.1f}")

    # The overall reaction is faster (total reactivity is higher)
    print(f"\nOverall reactivity increased from {total_reactivity:.1f} to {total_reactivity_excited:.1f}.")

    # The selectivity is dramatically changed.
    selectivity_ratio = prob_H_removal_excited / prob_D_removal_excited
    print(f"Likelihood of H removal is now {selectivity_ratio:.2f} times that of D removal.")

    print("\nConclusion: Exciting the C-H bond accelerates the reaction and vastly enhances the likelihood of H atom removal over D atoms.")

# Run the model
model_bond_selective_chemistry()