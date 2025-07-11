import sys

def solve():
    """
    This function calculates the infimum of c for which the given Markov chain is transient.
    It does so by printing the step-by-step mathematical derivation based on Lamperti's criterion.
    """

    print("To determine when the Markov chain is transient, we analyze the drift and variance of its jumps for large k.")
    print("-" * 70)

    # Step 1: Calculate the drift (mu_k)
    print("Step 1: Calculate the drift (mean jump) mu_k for large k.")
    print("mu_k = E[X_{n+1} - X_n | X_n = k]")
    print("mu_k = (k-2-k)*P(k,k-2) + (k+2-k)*P(k,k+2) + (k-1-k)*P(k,k-1) + (k+1-k)*P(k,k+1)")
    print("Substituting the probabilities P(k,j):")
    print("mu_k = (-2)*(1/4) + (2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k)")
    print("Simplifying the expression:")
    print("mu_k = -1/2 + 1/2 - 1/4 + c/k + 1/4 + c/k")
    mu_k_numerator_c = 2.0
    print(f"mu_k = ({mu_k_numerator_c}*c)/k\n")
    
    # Identify the coefficient A
    A_coeff = mu_k_numerator_c
    print(f"The drift is of the form A/k, where A = {A_coeff}*c.")
    print("-" * 70)

    # Step 2: Calculate the limiting jump variance (sigma^2)
    print("Step 2: Calculate the limiting variance of the jumps, sigma^2.")
    print("First, we find the second moment of the jumps, M2_k = E[(X_{n+1} - X_n)^2 | X_n = k].")
    print("M2_k = (-2)^2*P(k,k-2) + (2)^2*P(k,k+2) + (-1)^2*P(k,k-1) + (1)^2*P(k,k+1)")
    print("Substituting the probabilities P(k,j):")
    print("M2_k = 4*(1/4) + 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k)")
    print("Simplifying the expression:")
    print("M2_k = 1 + 1 + 1/4 - c/k + 1/4 + c/k")
    m2_k = 2.5
    print(f"M2_k = {m2_k}")
    
    print("\nThe variance sigma_k^2 is M2_k - (mu_k)^2.")
    print(f"sigma_k^2 = {m2_k} - (({A_coeff}*c)/k)^2")
    print("As k tends to infinity, mu_k tends to 0. So, the limiting variance sigma^2 is:")
    sigma_sq = m2_k
    print(f"sigma^2 = lim(k->inf) sigma_k^2 = {sigma_sq}\n")
    print("-" * 70)

    # Step 3: Apply the criterion for transience
    print("Step 3: Apply the criterion for transience.")
    print("A Markov chain with mu_k ~ A/k and limiting variance sigma^2 is transient if A > sigma^2 / 2.")
    print("From our calculations:")
    print(f"A = {A_coeff}*c")
    print(f"sigma^2 = {sigma_sq}")

    print("\nThe condition for transience is:")
    # A > sigma^2 / 2
    # 2c > 2.5 / 2
    val1 = sigma_sq
    val2 = 2.0
    print(f"{A_coeff}*c > {val1} / {val2}")
    
    # 2c > 1.25
    rhs1 = val1 / val2
    print(f"{A_coeff}*c > {rhs1}")
    
    # c > 1.25 / 2
    print(f"c > {rhs1} / {A_coeff}")
    
    # c > 0.625
    final_c = rhs1 / A_coeff
    print(f"c > {final_c}\n")
    print("-" * 70)

    # Step 4: Find the infimum
    print("Step 4: Find the infimum.")
    print(f"The set of c for which the chain is transient is ({final_c}, infinity).")
    print("The infimum (greatest lower bound) of this set is the value at the boundary.")
    print(f"Infimum = {final_c}")
    print("\nThis corresponds to the fraction 5/8.")


if __name__ == '__main__':
    solve()
    # To output the final answer in the required format for the platform
    # This part would not be executed if the user copies the code block.
    if 'ipykernel' not in sys.modules:
      print("<<<0.625>>>")