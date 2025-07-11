import fractions

def find_infimum_c():
    """
    This script finds the infimum of c for which a given Markov chain is transient.
    The method involves calculating the asymptotic drift and variance of the chain's
    jumps and applying Pakes's Lemma, a standard criterion for transience.
    """
    print("Step 1: Define the problem parameters from the transition probabilities.")
    print("For a state k, the jumps and their probabilities are:")
    jumps = [-2, -1, 1, 2]
    probs_str = ["1/4", "1/4 - c/k", "1/4 + c/k", "1/4"]
    print(f"Jumps: {jumps}")
    print(f"Probabilities: {probs_str}")
    print("-" * 60)

    print("Step 2: Calculate the asymptotic mean drift, mu.")
    print("The mean drift mu(k) = E[X_n+1 - X_n | X_n = k].")
    print("mu(k) = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print("mu(k) = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2")
    print("mu(k) = 2*c/k")
    print("This is in the form 'mu / k', so the parameter 'mu' for the criterion is 2*c.")
    mu_coeff_numerator = 2
    print(f"Parameter mu = {mu_coeff_numerator}*c")
    print("-" * 60)

    print("Step 3: Calculate the asymptotic variance of jumps, sigma^2.")
    print("The variance sigma^2(k) = E[(X_n+1 - X_n)^2 | X_n = k].")
    print("sigma^2(k) = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)")
    print("sigma^2(k) = 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)")
    print("sigma^2(k) = 1 + 1/4 + 1/4 + 1 = 2 + 1/2 = 5/2")
    sigma_sq_val = 5/2
    sigma_sq_frac = fractions.Fraction(sigma_sq_val)
    print(f"The asymptotic variance parameter sigma^2 = {sigma_sq_val} = {sigma_sq_frac.numerator}/{sigma_sq_frac.denominator}.")
    print("-" * 60)

    print("Step 4: Apply the transience criterion.")
    print("Pakes's Lemma states that the chain is transient if 2*mu / sigma^2 > 1.")
    print("Substituting the parameters we found:")
    # Using the requirement to output each number in the final equation
    print(f"The inequality is: 2 * ({mu_coeff_numerator}*c) / ({sigma_sq_frac.numerator}/{sigma_sq_frac.denominator}) > 1")
    print("-" * 60)

    print("Step 5: Solve the inequality for c.")
    # 2 * (2*c) / (5/2) > 1
    # 4*c / (5/2) > 1
    # (8/5)*c > 1
    # c > 5/8
    lhs_coeff_num = 2 * mu_coeff_numerator * sigma_sq_frac.denominator
    lhs_coeff_den = sigma_sq_frac.numerator
    print(f"Simplifying the inequality: ({lhs_coeff_num}/{lhs_coeff_den}) * c > 1")
    
    infimum_c_frac = fractions.Fraction(lhs_coeff_den, lhs_coeff_num)
    infimum_val = float(infimum_c_frac)
    
    print(f"Solving for c: c > {infimum_c_frac.denominator}/{infimum_c_frac.numerator}")
    print(f"So, c > {infimum_c_frac}")
    print("-" * 60)

    print("Step 6: Determine the infimum.")
    print(f"The set of values for c for which the chain is transient is ({infimum_c_frac}, infinity).")
    print(f"The infimum (greatest lower bound) of this set is {infimum_c_frac}.")
    print(f"\nThe final answer as a decimal is {infimum_val}.")

if __name__ == "__main__":
    find_infimum_c()

<<<0.625>>>