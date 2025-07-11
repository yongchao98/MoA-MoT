import math

def derive_photon_gas_properties():
    """
    This script outlines the derivation of equilibrium energy and entropy for a
    photon gas using the principles of statistical mechanics.
    """

    print("### Derivation of Equilibrium Properties for a Photon Gas (Bose Case) ###\n")

    # --- Step 1: The Principle of Maximum Entropy ---
    print("Step 1: The Principle of Maximum Entropy")
    print("----------------------------------------")
    print("The equilibrium state of an isolated system is the one that maximizes entropy (S).")
    print("This is the foundational concept, linking to large deviation theory: the equilibrium")
    print("state is the most probable one, and all other states are 'large deviations'.")
    print("For bosons, the entropy S is given by S = k_B * ln(Omega), where Omega is the")
    print("number of ways to distribute n_i particles into g_i states at energy level e_i:")
    print("S/k_B = sum_i [ (n_i + g_i)ln(n_i + g_i) - n_i*ln(n_i) - g_i*ln(g_i) ]\n")

    # --- Step 2: Constraints and the Bose-Einstein Distribution ---
    print("Step 2: Constraints and Equilibrium Distribution")
    print("-----------------------------------------------")
    print("We maximize S subject to one constraint for a photon gas: fixed total energy.")
    print("Constraint: E = sum_i(n_i * e_i)")
    print("(Note: Photon number N is not conserved, so its chemical potential is zero.)")
    print("Using Lagrange multipliers, this maximization yields the Bose-Einstein distribution for photons:")
    print("This gives the mean occupation number <n_i> per state (g_i=1) at energy e_i.")
    print("Let beta = 1 / (k_B * T).\n")

    print(">>> The Bose-Einstein distribution for photons is:")
    print("    <n(e)> = 1 / (exp(beta * e) - 1)")
    print("    <n(e)> = 1 / (exp(e / (k_B * T)) - 1)\n")


    # --- Step 3: Mean Energy Calculation ---
    print("Step 3: Calculating the Mean Energy (E)")
    print("---------------------------------------")
    print("The total energy E is the integral of energy 'e' times the number of photons")
    print("n(e) over all states. We use the density of states g(e) for photons in a volume V:")
    print("g(e) = (8 * pi * V * e^2) / (h^3 * c^3)\n")
    print("E = integral from 0 to infinity of [ e * <n(e)> * g(e) ] de")
    print("E = integral [ e * (1/(exp(beta*e) - 1)) * (8*pi*V*e^2 / (h^3*c^3)) ] de")
    print("This simplifies to an expression containing a standard integral:")
    print("Integral part: integral from 0 to inf [x^3 / (e^x - 1)] dx = pi^4 / 15\n")

    print(">>> The final equilibrium Mean Energy E is given by the Stefan-Boltzmann Law:")
    energy_numerator = "8 * pi^5 * V * k_B^4 * T^4"
    energy_denominator = "15 * h^3 * c^3"
    print(f"    E = ({energy_numerator}) / ({energy_denominator})\n")
    print("Breaking down the numbers in the final equation for Energy:")
    print(f"    Numerical constant in numerator: 8")
    print(f"    Power of pi: 5")
    print(f"    Power of Boltzmann constant k_B: 4")
    print(f"    Power of Temperature T: 4")
    print(f"    Numerical constant in denominator: 15")
    print(f"    Power of Planck constant h: 3")
    print(f"    Power of speed of light c: 3\n")


    # --- Step 4: Entropy Calculation ---
    print("Step 4: Calculating the Equilibrium Entropy (S)")
    print("-----------------------------------------------")
    print("Entropy can be found from the thermodynamic relation S = (E - F) / T,")
    print("where F is the Helmholtz Free Energy.")
    print("For a photon gas, it can be shown that F = -E / 3.\n")
    print("Therefore, S = (E - (-E/3)) / T = (4/3) * E / T\n")
    print(">>> The final equilibrium Entropy S is:")

    entropy_numerator = f"4 * ({energy_numerator})"
    entropy_denominator = f"3 * T * ({energy_denominator})"
    final_entropy_numerator = "32 * pi^5 * V * k_B^4 * T^3"
    final_entropy_denominator = "45 * h^3 * c^3"

    print(f"    S = (4/3) * E / T\n")
    print("Substituting the expression for E, we get:")
    print(f"    S = ({final_entropy_numerator}) / ({final_entropy_denominator})\n")
    print("Breaking down the numbers in the final equation for Entropy:")
    print(f"    Numerical constant in numerator: 32 (from 4*8)")
    print(f"    Power of pi: 5")
    print(f"    Power of Boltzmann constant k_B: 4")
    print(f"    Power of Temperature T: 3")
    print(f"    Numerical constant in denominator: 45 (from 3*15)")
    print(f"    Power of Planck constant h: 3")
    print(f"    Power of speed of light c: 3\n")

if __name__ == "__main__":
    derive_photon_gas_properties()