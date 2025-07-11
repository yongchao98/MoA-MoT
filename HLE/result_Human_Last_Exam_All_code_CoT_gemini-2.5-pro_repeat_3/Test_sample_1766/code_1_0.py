from fractions import Fraction

def solve_spin_orbit_coupling():
    """
    Calculates and prints the J^2 and Jz common eigenvalues for a p-electron.
    """
    # Step 1: Define quantum numbers for a p-electron
    l = 1  # orbital angular momentum for p-shell
    s = 0.5  # spin for an electron

    print(f"For a p-electron, the quantum numbers are l = {l} and s = {Fraction(s)}.")

    # Step 2: Determine possible j values
    j_min = abs(l - s)
    j_max = l + s
    
    # In this case, j can be |1-1/2|=1/2 and 1+1/2=3/2
    j_values = [j_min, j_max] 
    
    j_values_str = ', '.join([str(Fraction(j)) for j in j_values])
    print(f"The possible total angular momentum quantum numbers are j = {j_values_str}.\n")

    # Step 3: Calculate eigenvalues for each j
    for j in j_values:
        j_frac = Fraction(j)
        j2_eigenvalue = j * (j + 1)
        j2_eigenvalue_frac = Fraction(j2_eigenvalue)
        
        print(f"--- For the case j = {j_frac} ---")
        
        # Calculate possible m_j values
        m_j_values = [x * 0.5 for x in range(int(-j * 2), int(j * 2) + 2, 2)]
        mj_values_str = ', '.join([str(Fraction(mj)) for mj in m_j_values])
        print(f"The possible magnetic quantum numbers m_j are: {mj_values_str}.")
        
        print("\nThe eigenvalue equations for the common eigenstates |j, m_j> are:")
        
        # Print the final equations for each state
        for m_j in m_j_values:
            m_j_frac = Fraction(m_j)
            
            # Print state information
            print(f"For the state |{j_frac}, {m_j_frac}>:")
            
            # Print J^2 eigenvalue equation
            print(f"  J²|{j_frac}, {m_j_frac}> = {j2_eigenvalue_frac}ħ² |{j_frac}, {m_j_frac}>")
            
            # Print J_z eigenvalue equation
            print(f"  Jz|{j_frac}, {m_j_frac}> = {m_j_frac}ħ |{j_frac}, {m_j_frac}>\n")

solve_spin_orbit_coupling()
<<<The possible j values are 1/2 and 3/2.>>>