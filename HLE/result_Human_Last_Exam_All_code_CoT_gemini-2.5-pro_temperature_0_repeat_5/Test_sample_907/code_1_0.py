import numpy as np

def print_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules under two different assumptions.
    """

    print("----------------------------------------------------------------------")
    print("Equation for Absorption Cross-Section (First-Order Perturbation Theory)")
    print("----------------------------------------------------------------------")
    print("The general form of the absorption cross-section σ for a system interacting")
    print("with a Gaussian laser pulse of frequency ω_L and duration τ is:")
    print("\n  σ(ω_L) = C * Σ_f |<f|D|i>|² * ω_fi * exp(-(ω_fi - ω_L)² * τ²)\n")
    print("Where:")
    print("  C     = A proportionality constant")
    print("  |i>   = Initial ground state of the molecular chain")
    print("  |f>   = Final exciton state")
    print("  D     = Total electric dipole operator of the chain")
    print("  ω_fi  = Transition frequency from state |i> to |f>")
    print("  τ     = Duration of the Gaussian laser pulse")
    print("  ω_L   = Central frequency of the laser pulse")
    print("----------------------------------------------------------------------\n")

    # --- Case a) No interaction between molecules ---
    print("a) Case: The interaction between molecules can be neglected.")
    print("----------------------------------------------------------------------")
    print("In this case, all N molecules are independent and identical. The excitation")
    print("is localized on a single molecule. There is only one transition frequency ω_0.")
    print("The total cross-section is N times the cross-section of a single molecule.")
    print("\nThe equation is:")
    print("\n  σ(ω_L) = C * N * |μ|² * ω_0 * exp(-(ω_0 - ω_L)² * τ²)\n")
    print("Where:")
    print("  N     = Number of molecules in the chain")
    print("  μ     = Transition dipole moment of a single molecule")
    print("  ω_0   = Transition frequency of a single molecule")
    print("----------------------------------------------------------------------\n")

    # --- Case b) Near-neighbor interaction considered ---
    print("b) Case: The interaction between near-neighbors should be considered.")
    print("----------------------------------------------------------------------")
    print("Interaction between molecules leads to the formation of delocalized Frenkel")
    print("exciton states |Ψ_k>, where k = 1, 2, ..., N. The absorption spectrum")
    print("is a sum over transitions to these exciton states.")
    print("\nThe general equation is:")
    print("\n  σ(ω_L) = C * Σ_{k=1 to N} |d_k|² * ω_k * exp(-(ω_k - ω_L)² * τ²)\n")
    print("Where the sum is over the N exciton states, and:")
    print("\n  1. The transition frequency ω_k for each exciton state |Ψ_k> is:")
    print("     ω_k = ω_0 + (2J/ħ) * cos(π*k / (N+1))")
    print("     (J is the near-neighbor coupling energy, ħ is the reduced Planck constant)")
    print("\n  2. The transition dipole moment d_k for each exciton state is:")
    print("     d_k = <Ψ_k|D|i> = μ * sqrt(2/(N+1)) * Σ_{n=1 to N} sin(π*k*n / (N+1))")
    print("     (A selection rule applies: d_k is non-zero only for odd values of k)")
    print("\nLet's expand the equation for a chain of N=4 molecules to see its structure.")
    print("For N=4, the allowed transitions are to exciton states k=1 and k=3.")
    
    # Symbolic expansion for N=4
    N = 4
    # k=1 term
    omega_1_str = "ω_0 + (2J/ħ) * cos(π/5)"
    # Sum sin(pi*n/5) for n=1,2,3,4
    sum_sin_1 = np.sin(np.pi/5) + np.sin(2*np.pi/5) + np.sin(3*np.pi/5) + np.sin(4*np.pi/5)
    d_1_sq_str = f"|μ|² * (2/5) * ({sum_sin_1:.3f})²"

    # k=3 term
    omega_3_str = "ω_0 + (2J/ħ) * cos(3π/5)"
    # Sum sin(3*pi*n/5) for n=1,2,3,4
    sum_sin_3 = np.sin(3*np.pi/5) + np.sin(6*np.pi/5) + np.sin(9*np.pi/5) + np.sin(12*np.pi/5)
    d_3_sq_str = f"|μ|² * (2/5) * ({sum_sin_3:.3f})²"

    print("\nFor N=4, the final equation has two terms (k=1 and k=3):")
    print("\n  σ(ω_L) = C * [ Term(k=1) + Term(k=3) ]")
    
    print("\n  Term(k=1) = |d₁|² * ω₁ * exp(-(ω₁ - ω_L)² * τ²)")
    print(f"    where ω₁ = {omega_1_str}")
    print(f"          |d₁|² = {d_1_sq_str}")

    print("\n  Term(k=3) = |d₃|² * ω₃ * exp(-(ω₃ - ω_L)² * τ²)")
    print(f"    where ω₃ = {omega_3_str}")
    print(f"          |d₃|² = {d_3_sq_str}")
    print("----------------------------------------------------------------------")

if __name__ == '__main__':
    print_absorption_equations()
    # The final answer is the explanation and set of equations derived.
    # The python code above generates this explanation.
    # As no single value is requested, the full textual output is the answer.
    # To comply with the format, I will output a placeholder.
    final_answer = "<<<The final answer is the complete derivation and set of equations provided by the Python script's output.>>>"
