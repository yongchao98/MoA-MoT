import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture of acids and a base.
    
    The problem involves:
    - 500 mL of 0.1 M CH3COOH (weak acid)
    - 400 mL of 0.2 M HCl (strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (strong base)
    
    The provided answer is 12.62, which corresponds to option A.
    """
    
    # Define initial conditions
    vol_ch3cooh = 0.500  # L
    mol_ch3cooh = 0.1    # M

    vol_hcl = 0.400      # L
    mol_hcl = 0.2        # M

    vol_baoh2 = 0.300    # L
    mol_baoh2 = 0.3      # M

    # Step 1: Calculate the initial moles of all acidic and basic species.
    moles_ch3cooh = vol_ch3cooh * mol_ch3cooh
    moles_h_from_hcl = vol_hcl * mol_hcl
    
    # Ba(OH)2 is a strong base that provides 2 OH- ions per formula unit.
    moles_oh_from_baoh2 = vol_baoh2 * mol_baoh2 * 2

    # Step 2: Determine the total moles of acid and base to be neutralized.
    # The strong base will react with all available protons, from both strong and weak acids.
    total_moles_acid = moles_h_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_from_baoh2

    # Step 3: Determine the excess reactant.
    if total_moles_base > total_moles_acid:
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # Step 4: Calculate the total volume of the solution.
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        
        # Step 5: Calculate the final concentration of the excess species.
        final_oh_concentration = excess_moles_oh / total_volume
        
        # Step 6: Calculate pOH and then pH.
        try:
            poh = -math.log10(final_oh_concentration)
            ph = 14 - poh
        except ValueError:
            return "Error: Cannot calculate log of a non-positive concentration."

    elif total_moles_acid > total_moles_base:
        # This case is more complex as it would leave a mixture of a strong acid and a weak acid.
        # However, based on the numbers, we expect the base to be in excess.
        # For completeness, let's calculate it.
        # The strong base would first neutralize the strong acid.
        remaining_h_from_hcl = moles_h_from_hcl - total_moles_base
        if remaining_h_from_hcl < 0: # Strong acid is fully neutralized
            # This path shouldn't be taken based on our initial check, but as a safeguard:
            return "Incorrect logic path: Base should be in excess."
        
        # This case is not applicable here, but the pH would be dominated by the remaining strong acid.
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_h_concentration = remaining_h_from_hcl / total_volume
        try:
            ph = -math.log10(final_h_concentration)
        except ValueError:
            return "Error: Cannot calculate log of a non-positive concentration."
            
    else: # Moles of acid and base are equal
        # This would result in a solution of the conjugate base of the weak acid (CH3COO-).
        # The pH would be basic, but this is not the case here.
        ph = 7.0 # Neutralization of strong acid/base, but weak acid's conjugate base remains.
                 # This is an oversimplification and not the correct path for this problem.

    # The expected answer is 12.62 (Option A)
    expected_ph = 12.62
    
    # Check if the calculated pH matches the expected answer within a small tolerance.
    if abs(ph - expected_ph) < 0.01:
        return "Correct"
    else:
        return (f"Incorrect. The calculated pH is {ph:.2f}, but the expected answer is {expected_ph}. "
                f"The calculation steps are as follows:\n"
                f"1. Moles of H+ (total) = {total_moles_acid:.3f} mol\n"
                f"2. Moles of OH- = {total_moles_base:.3f} mol\n"
                f"3. Excess moles of OH- = {excess_moles_oh:.3f} mol\n"
                f"4. Total Volume = {total_volume:.1f} L\n"
                f"5. Final [OH-] = {final_oh_concentration:.5f} M\n"
                f"6. pOH = {-math.log10(final_oh_concentration):.2f}\n"
                f"7. pH = 14 - pOH = {ph:.2f}")

# Run the check
result = check_ph_calculation()
print(result)