import math

def solve_steady_state_probability():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process and prints the derivation steps.
    """
    
    # Define symbolic variables for clarity in the output
    lambda_sym = "λ"
    mu_sym = "μ"
    rho_sym = "ρ"
    pi_0_sym = "π_0"
    
    print("Step 1: Identify the birth and death rates.")
    print(f"The birth rate from state i to i+1 is λ_i = {lambda_sym} / (i + 1).")
    print(f"The death rate from state i to i-1 is μ_i = {mu_sym}.")
    print("-" * 50)
    
    print("Step 2: Use the detailed balance equations.")
    print("The detailed balance equation is: π_i * λ_i = π_{i+1} * μ_{i+1}.")
    print("This allows us to express π_{i+1} in terms of π_i:")
    print("π_{i+1} = π_i * (λ_i / μ_{i+1})")
    print("-" * 50)

    print("Step 3: Express each π_n in terms of π_0.")
    print("By applying the relation iteratively, we get a general formula:")
    print("π_n = π_0 * Π_{i=0}^{n-1} (λ_i / μ_{i+1})")
    print("\nSubstituting the given rates:")
    print(f"λ_i = {lambda_sym} / (i + 1)")
    print(f"μ_{i+1} = {mu_sym}")
    print(f"π_n = {pi_0_sym} * Π_{{i=0}}^{{n-1}} (({lambda_sym} / (i + 1)) / {mu_sym})")
    print(f"π_n = {pi_0_sym} * Π_{{i=0}}^{{n-1}} ({lambda_sym} / ({mu_sym} * (i + 1)))")
    
    print(f"\nLet {rho_sym} = {lambda_sym} / {mu_sym}. The expression simplifies to:")
    print(f"π_n = {pi_0_sym} * Π_{{i=0}}^{{n-1}} ({rho_sym} / (i + 1))")
    print("Expanding the product term (Π):")
    print(f"π_n = {pi_0_sym} * ({rho_sym}/1) * ({rho_sym}/2) * ({rho_sym}/3) * ... * ({rho_sym}/n)")
    print(f"This simplifies to π_n = {pi_0_sym} * ({rho_sym}^n / n!)")
    print("-" * 50)
    
    print("Step 4: Apply the normalization condition Σ π_n = 1.")
    print(f"The sum of all probabilities must be 1: Σ_{{n=0}}^{{∞}} {pi_0_sym} * ({rho_sym}^n / n!) = 1")
    print(f"{pi_0_sym} * Σ_{{n=0}}^{{∞}} ({rho_sym}^n / n!) = 1")
    print("-" * 50)

    print("Step 5: Recognize the series and solve for π_0.")
    print("The series Σ_{n=0}^{∞} (x^n / n!) is the Taylor series expansion for e^x.")
    print(f"Therefore, the sum Σ_{{n=0}}^{{∞}} ({rho_sym}^n / n!) is equal to e^{rho_sym}.")
    print("\nSubstituting this into our equation:")
    print(f"{pi_0_sym} * e^{rho_sym} = 1")
    
    print("\nSolving for π_0 gives the final expression:")
    print(f"{pi_0_sym} = 1 / e^{rho_sym}")
    print("Which can be written as:")
    
    # Final equation output
    final_base = "e"
    final_exponent_coeff = -1
    final_exponent_var = "ρ"
    
    print(f"{pi_0_sym} = {final_base}^({final_exponent_coeff}*{final_exponent_var})")
    print(f"Or more simply: {pi_0_sym} = e^(-{rho_sym})")

if __name__ == "__main__":
    solve_steady_state_probability()