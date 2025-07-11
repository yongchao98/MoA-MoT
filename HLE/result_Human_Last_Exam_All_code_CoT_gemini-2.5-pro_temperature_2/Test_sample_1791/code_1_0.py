import sys

def derive_effective_interaction():
    """
    This script outlines the derivation of the effective electron-electron 
    interaction by integrating out phonon modes and prints the final formula.
    """

    # For display purposes, use unicode for nicer symbols where possible
    hbar = "\u210F"
    sum_symbol = "\u2211"
    rho = "\u03C1"
    omega = "\u03C9"
    nu = "\u03BD"
    delta = "\u03B4"
    tau = "\u03C4"
    
    if sys.stdout.encoding != 'utf-8':
        # Fallback to ASCII if shell does not support UTF-8
        hbar = "w" # user wrote w_q for frequency, probably meant omega and set hbar=1
        sum_symbol = "sum"
        rho = "rho"
        omega = "w"
        nu = "nu"
        delta = "delta"
        tau = "tau"


    print("### Step 1: Defining the System ###")
    print("-----------------------------------")
    print(f"The phonon Hamiltonian is given by H_ph = {sum_symbol}_{{q,j}} {omega}_q a_{{q,j}}\u207A a_{{q,j}}.")
    print(f"The electron density operator is {rho}(q).")
    print("The given electron-phonon coupling is:")
    print(f"H_el-ph = {sum_symbol}_{{q,j}} C_{{q,j}} {rho}_{{q}} (a_{{q,j}} + a_{{-q,j}}\u207A),")
    print(f"with the coupling constant C_{{q,j}} = g * i * q_j / (2 * m * {omega}_q)^(1/2).\n")

    print("### Step 2: Integrating out Phonons ###")
    print("-----------------------------------------")
    print("To find the effective interaction between electrons, we integrate out the phonon fields.")
    print("This is equivalent to calculating the change in the electronic action to second order in H_el-ph:")
    print(f"Delta S_eff = -1/2 * \u222B d{tau}1 d{tau}2 <T_{tau} H_el-ph({tau}1) H_el-ph({tau}2)>_ph")
    print("where <...>_ph denotes the thermal average over the free phonon system.\n")
    
    print("### Step 3: Calculating the Phonon Propagator ###")
    print("-------------------------------------------------")
    print(f"The expectation value requires evaluating the phonon Green's function. The interaction term contains operators of the form X_q = (a_q + a_{-q}\u207A).")
    print(f"We need to compute the propagator D_ph({tau}1-{tau}2) = <T_{tau} X_{{q,j}}({tau}1) X_{{-q,j}}({tau}2)>_ph.")
    print(f"This is because modes are independent (j'=j) and momentum must be conserved (q'=-q).")
    print("In Matsubara frequency space, this propagator evaluates to:")
    print(f"D_ph(q, j, i{nu}_n) = 2 * {omega}_q / ((i{nu}_n)^2 - {omega}_q^2) = -2 * {omega}_q / ({nu}_n^2 + {omega}_q^2).\n")

    print("### Step 4: Deriving the Effective Potential ###")
    print("------------------------------------------------")
    print("The effective action takes the form:")
    print(f"Delta S_eff = 1/2 * (1/beta) * {sum_symbol}_{{q,n}} U_eff(q, i{nu}_n) {rho}_q(i{nu}_n) {rho}_{-q}(-i{nu}_n).")
    print("By evaluating the second-order expression, we get:")
    print(f"Delta S_eff = -1/2 * (1/beta) * {sum_symbol}_{{q,j,n}} [C_{{q,j}} C_{{-q,j}}] D_ph(q, j, i{nu}_n) {rho}_q(i{nu}_n) {rho}_{-q}(-i{nu}_n).")
    print("The product of the coupling constants is:")
    print( "C_{q,j} * C_{-q,j} = [g*i*q_j / sqrt(2m{omega}_q)] * [g*i*(-q_j) / sqrt(2m{omega}_q)]")
    print(f"                 = g^2 * (i^2) * (-q_j^2) / (2*m*{omega}_q) = g^2 * q_j^2 / (2*m*{omega}_q).")
    print("\nThus, the contribution from polarization 'j' to the potential is:")
    print(f"U_j(q, i{nu}_n) = -[C_{{q,j}} C_{{-q,j}}] * D_ph(q, j, i{nu}_n)")
    print(f"             = -[g^2 * q_j^2 / (2*m*{omega}_q)] * [-2*{omega}_q / ({nu}_n^2 + {omega}_q^2)]")
    print(f"             = g^2 * q_j^2 / (m * ({omega}_q^2 + {nu}_n^2)).")
    print("\nThe total potential is the sum over polarizations j (e.g., j=x,y,z):")
    print(f"U_eff(q, i{nu}_n) = {sum_symbol}_j U_j = {sum_symbol}_j g^2 * q_j^2 / (m * ({omega}_q^2 + {nu}_n^2))")
    print(f"Assuming {sum_symbol}_j q_j^2 = q_x^2 + q_y^2 + q_z^2 = |q|^2, we get the final result.\n")

    print("### Step 5: The Final Equation ###")
    print("----------------------------------")
    print("The effective electron-electron interaction potential U_eff(q) is frequency-dependent:")

    g = "g"
    q_vec_sq = "q^2"
    m = "m"
    w_q_sq = f"{omega}_q^2"
    nu_n_sq = f"{nu}_n^2"

    numerator_str = f"g^2 * q^2"
    denominator_str = f"m * ( {w_q_sq} + {nu_n_sq} )"

    print(f"U_eff(q, i{nu}_n) = {numerator_str} / ({denominator_str})")
    
    print("\nIndividual parts of the equation:")
    print(f"g^2:       The square of the electron-phonon coupling constant.")
    print(f"q^2:       The square of the magnitude of the momentum transfer vector q.")
    print(f"m:         The mass parameter from the coupling definition.")
    print(f"{omega}_q^2:     The square of the phonon frequency at momentum q.")
    print(f"{nu}_n^2:     The square of the bosonic Matsubara frequency.")
    
    
derive_effective_interaction()
<<<U_eff(q, i*nu_n) = g^2 * q^2 / (m * (w_q^2 + nu_n^2))>>>