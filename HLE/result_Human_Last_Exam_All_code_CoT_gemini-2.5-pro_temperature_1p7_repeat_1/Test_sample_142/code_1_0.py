import math

def calculate_csf():
    """
    Calculates the number of singlet CSFs for a Full CI calculation
    of CH2SiHF with a 6-31G** basis set.
    """
    # Step 1: Define the number of electrons.
    # C(6) + H(3) + Si(14) + F(9) = 32 electrons.
    n_elec = 32
    # For a singlet state, n_alpha = n_beta = n_elec / 2
    n_alpha = n_elec // 2

    # Step 2: Define the number of basis functions (orbitals).
    # C, Si, F in 6-31G** are 15 basis functions each.
    # H in 6-31G** is 5 basis functions each.
    # Total = 1*C + 3*H + 1*Si + 1*F = 1*15 + 3*5 + 1*15 + 1*15 = 60
    n_orb = 60

    # Step 3: Apply the formula for the number of singlet CSFs.
    # Formula: N_csf = (comb(n_orb, n_alpha))^2 - (comb(n_orb, n_alpha - 1))^2
    
    # Calculate the binomial coefficients
    term1_comb = math.comb(n_orb, n_alpha)
    term2_comb = math.comb(n_orb, n_alpha - 1)
    
    # Calculate the number of CSFs
    n_csf = term1_comb**2 - term2_comb**2

    # Print the equation with the calculated numbers
    print(f"The number of orbitals (N_orb) is: {n_orb}")
    print(f"The number of alpha electrons (N_alpha) is: {n_alpha}")
    print("\nCalculating the number of singlet CSFs using the formula:")
    print(f"N_CSF = (C(N_orb, N_alpha))^2 - (C(N_orb, N_alpha - 1))^2")
    print(f"N_CSF = (C({n_orb}, {n_alpha}))^2 - (C({n_orb}, {n_alpha - 1}))^2")
    print(f"N_CSF = ({term1_comb})^2 - ({term2_comb})^2")
    print(f"\nThe total number of singlet CSFs is: {n_csf}")

calculate_csf()