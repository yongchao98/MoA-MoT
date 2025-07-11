def solve_steady_state_probability():
    """
    This script explains the derivation of the steady-state probability pi_0
    for the given birth-death process.
    """
    # Define symbolic variables for clarity in the explanation.
    rho_str = "ρ"
    lambda_str = "λ"
    mu_str = "μ"
    pi_0_str = "π_0"
    pi_i_str = "π_i"
    pi_i_plus_1_str = "π_{i+1}"
    pi_n_str = "π_n"
    e_str = "e"

    print("Step 1: Set up the detailed balance equations.")
    print(f"The general equation for a birth-death process is: {pi_i_str} * λ_i = {pi_i_plus_1_str} * μ_{i+1}")
    print(f"For this process, the birth rate is λ_i = {lambda_str} / (i+1) and the death rate is μ_{i+1} = {mu_str}.")
    print(f"Substituting these rates, we get: {pi_i_str} * ({lambda_str} / (i+1)) = {pi_i_plus_1_str} * {mu_str}")
    print("-" * 30)

    print("Step 2: Express π_{i+1} in terms of π_i.")
    print(f"Rearranging the equation: {pi_i_plus_1_str} = {pi_i_str} * ({lambda_str} / ({mu_str} * (i+1)))")
    print(f"Let {rho_str} = {lambda_str} / {mu_str}. The equation becomes: {pi_i_plus_1_str} = {pi_i_str} * ({rho_str} / (i+1))")
    print("-" * 30)

    print("Step 3: Express all π_n in terms of π_0 by induction.")
    print(f"For n=1: π_1 = π_0 * ({rho_str} / 1) = {pi_0_str} * {rho_str}")
    print(f"For n=2: π_2 = π_1 * ({rho_str} / 2) = ({pi_0_str} * {rho_str}) * ({rho_str} / 2) = {pi_0_str} * ({rho_str}^2 / 2!)")
    print(f"For n=3: π_3 = π_2 * ({rho_str} / 3) = ({pi_0_str} * {rho_str}^2 / 2!) * ({rho_str} / 3) = {pi_0_str} * ({rho_str}^3 / 3!)")
    print(f"The general formula emerges: {pi_n_str} = {pi_0_str} * ({rho_str}^n / n!)")
    print("-" * 30)

    print("Step 4: Use the normalization condition.")
    print("The sum of all probabilities must equal 1: Σ (from n=0 to ∞) π_n = 1")
    print(f"Substituting our formula for {pi_n_str}: Σ (from n=0 to ∞) {pi_0_str} * ({rho_str}^n / n!) = 1")
    print(f"Factoring out {pi_0_str}: {pi_0_str} * Σ (from n=0 to ∞) ({rho_str}^n / n!) = 1")
    print("-" * 30)

    print("Step 5: Recognize the mathematical series.")
    print(f"The summation Σ ({rho_str}^n / n!) is the Taylor series expansion for the exponential function e^x, where x = {rho_str}.")
    print(f"So, the summation equals {e_str}^{rho_str}.")
    print(f"The normalization equation simplifies to: {pi_0_str} * {e_str}^{rho_str} = 1")
    print("-" * 30)

    print("Step 6: Solve for π_0.")
    print("The final expression for the steady-state probability π_0 is:")
    final_equation = f"{pi_0_str} = 1 / ({e_str}^{rho_str}) = {e_str}^(-{rho_str})"
    print(final_equation)

if __name__ == '__main__':
    solve_steady_state_probability()