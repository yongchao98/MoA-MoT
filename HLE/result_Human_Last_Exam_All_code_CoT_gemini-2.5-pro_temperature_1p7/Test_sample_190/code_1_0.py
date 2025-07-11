from fractions import Fraction

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for a Markov chain to be transient.

    The method uses Lamperti's criterion for random walks on Z+.
    A chain is transient if alpha > B/2, where:
    - alpha = lim_{k->inf} k * mu_k
    - B = lim_{k->inf} E[(X_{n+1} - X_{n})^2 | X_n = k]
    - mu_k = E[X_{n+1} - X_{n} | X_n = k]
    """

    print("Step 1: Define transition probabilities and jumps for large k.")
    # Transition probabilities:
    # P_{k,k-2} = 1/4
    # P_{k,k-1} = 1/4 - c/k
    # P_{k,k+1} = 1/4 + c/k
    # P_{k,k+2} = 1/4
    jumps_and_probs = {
        -2: "1/4",
        -1: "1/4 - c/k",
        1:  "1/4 + c/k",
        2:  "1/4"
    }
    print("Jumps and their probabilities:", jumps_and_probs)
    print("-" * 30)

    print("Step 2: Calculate the mean drift mu_k.")
    # mu_k = sum_{j} (j-k) * P_{k,j}
    # mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)
    # mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2
    # mu_k = 2*c/k
    print("The mean drift mu_k = (-2)(1/4) + (-1)(1/4 - c/k) + (1)(1/4 + c/k) + (2)(1/4) = 2*c/k.")
    # The calculation shows mu_k is proportional to 2c/k. The coefficient of c/k is 2.
    mu_k_coeff = 2
    print("-" * 30)
    
    print("Step 3: Calculate the parameter alpha.")
    # alpha = lim_{k->inf} k * mu_k = lim_{k->inf} k * (2*c/k) = 2c
    alpha_coeff = mu_k_coeff
    print(f"alpha = lim (k -> inf) k * mu_k = {alpha_coeff}*c")
    print("-" * 30)

    print("Step 4: Calculate the mean squared jump B.")
    # B = lim_{k->inf} E[(X_{n+1} - X_{n})^2 | X_n = k]
    # E_k = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)
    # E_k = 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)
    # E_k = 1 + 1/4 - c/k + 1/4 + c/k + 1 = 2.5
    B = 4*(1/4) + 1*(1/4) + 1*(1/4) + 4*(1/4)
    print("The mean squared jump B_k = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4) = 2.5")
    print(f"B = lim (k -> inf) B_k = {B}")
    print("-" * 30)

    print("Step 5: Apply the transience criterion: alpha > B / 2.")
    # The final equation is alpha_coeff * c > B / 2
    print("The inequality for transience is:")
    # Per instructions, printing each number in the final equation.
    print(f"{alpha_coeff} * c > {B} / 2")
    print("-" * 30)

    print("Step 6: Solve for c.")
    c_infimum = (B / 2) / alpha_coeff
    print(f"c > {B / 2} / {alpha_coeff}")
    print(f"c > {c_infimum}")
    print("-" * 30)

    print("Step 7: Find the infimum.")
    print(f"The set of c for which the chain is transient is ({c_infimum}, infinity).")
    print(f"The infimum (greatest lower bound) of this set is {c_infimum}.")
    
    # Final answer as a fraction
    infimum_fraction = Fraction(c_infimum).limit_denominator()
    print(f"The infimum in fractional form is {infimum_fraction}.")


if __name__ == "__main__":
    solve_markov_chain_transience()
    
    # Calculating the value to be put in the <<<>>> format
    alpha_coeff_val = 2
    B_val = 2.5
    c_infimum_val = (B_val / 2) / alpha_coeff_val
    final_answer = Fraction(c_infimum_val).limit_denominator()
    # This value will not be printed directly but will be used for the final submission format.
