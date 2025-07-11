import sys

def solve_partition_function():
    """
    This function explains and prints the formula for the grand canonical partition function
    Z using the path integral formalism.
    """
    # Use UTF-8 encoding for special characters if possible
    if sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
        hbar = 'h_bar'
        d_tau = 'd_tau'
        del_tau = 'd/d_tau'
        psi_star = 'psi*'
        psi = 'psi'
        mu = 'mu'
        beta = 'beta'
        integral_sign = 'INT'
    else:
        hbar = 'ħ'
        d_tau = 'dτ'
        del_tau = '∂_τ'
        psi_star = 'ψ*'
        psi = 'ψ'
        mu = 'μ'
        beta = 'β'
        integral_sign = '∫'

    print("Derivation of the Grand Canonical Partition Function Z via Path Integrals\n")

    print("Step 1: The Definition in Operator Formalism")
    print("---------------------------------------------")
    print(f"The grand canonical partition function (Z) describes a system at constant temperature (T) and chemical potential ({mu}). It is defined as the trace (Tr) over the system's states:")
    print(f"  Z = Tr[exp(-{beta} * (H - {mu}N))]")
    print(f"Where:")
    print(f"  H    : The Hamiltonian operator of the system.")
    print(f"  N    : The particle number operator.")
    print(f"  {mu}    : The chemical potential.")
    print(f"  {beta}   : Inverse temperature, 1 / (k_B * T).\n")

    print("Step 2: The Path Integral Representation")
    print("------------------------------------------")
    print("The trace can be evaluated by summing over all possible field configurations in imaginary time. This is the path integral.")
    print(f"  Z = {integral_sign} D[{psi_star}, {psi}] exp(-S_E[{psi_star}, {psi}])")
    print("Where:")
    print(f"  {integral_sign} D[{psi_star}, {psi}] : The functional integral over all field configurations.")
    print(f"  S_E         : The Euclidean action, which determines the 'weight' of each configuration.\n")

    print("Step 3: Finding the Euclidean Action S_E")
    print("------------------------------------------")
    print(f"The Euclidean action is found by taking the system's Lagrangian, including the chemical potential term, and transforming it from real time (t) to imaginary time ({d_tau}).")
    print(f"The term '-{mu}N' in the Hamiltonian effectively adds a '({mu}{psi_star}{psi})' term to the Lagrangian density.")
    print(f"This modifies the part of the action containing the time derivative. After the transformation (Wick Rotation), the chemical potential {mu} appears alongside the imaginary time derivative ({del_tau}).\n")

    print("Step 4: The Final Formula")
    print("--------------------------")
    print("The grand canonical partition function Z is the path integral over the fields with the following Euclidean action S_E:")

    # Printing the final equation with each "number" (component)
    print("\n--- FINAL EQUATION ---\n")
    final_equation = (f"  S_E[{psi_star}, {psi}] = {integral_sign}_0^{beta} {d_tau} {integral_sign} d^d x "
                      f"[ {psi_star}(x, {d_tau}) * ({hbar}{del_tau} - {mu}) * {psi}(x, {d_tau}) + H({psi_star}, {psi}) ]")
    print(final_equation)
    print("\n--- Equation Components ---")
    print(f"S_E       : The Euclidean Action.")
    print(f"{integral_sign}_0^{beta} {d_tau}  : Integral over imaginary time from 0 to {beta}.")
    print(f"{integral_sign} d^d x    : Integral over d spatial dimensions.")
    print(f"{psi_star}, {psi}    : The system's fields (e.g., complex numbers for bosons, Grassmann numbers for fermions).")
    print(f"{hbar}{del_tau}       : The kinetic term in imaginary time.")
    print(f"-{mu}         : The chemical potential, which shifts the energy levels.")
    print(f"H({psi_star}, {psi}) : The original Hamiltonian density (e.g., spatial kinetic energy and interaction terms).\n")

if __name__ == '__main__':
    solve_partition_function()