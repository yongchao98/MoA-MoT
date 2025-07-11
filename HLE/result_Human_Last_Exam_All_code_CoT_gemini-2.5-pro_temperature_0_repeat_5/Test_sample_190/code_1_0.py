import fractions

def find_infimum_c():
    """
    This script calculates the infimum of the set of values 'c' for which a
    given Markov chain is transient.

    The method involves analyzing the first and second moments of the chain's jumps.
    """
    print("Analyzing the Markov chain to find the condition for transience.")
    print("The transition probabilities for a large state k are given as:")
    print("P(k, k-2) = 1/4")
    print("P(k, k+2) = 1/4")
    print("P(k, k-1) = 1/4 - c/k")
    print("P(k, k+1) = 1/4 + c/k")
    print("-" * 50)

    # Step 1: Calculate the drift (mean jump) mu_k
    print("Step 1: Calculate the drift mu_k = E[X_{n+1} - X_n | X_n = k].")
    print("mu_k = (-2) * P(k, k-2) + (-1) * P(k, k-1) + (1) * P(k, k+1) + (2) * P(k, k+2)")
    print("mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print("mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2")
    print("After simplification, the drift is: mu_k = 2*c/k")
    print("-" * 50)

    # Step 2: Calculate the expected squared jump sigma_k^2
    print("Step 2: Calculate the expected squared jump sigma_k^2 = E[(X_{n+1} - X_n)^2 | X_n = k].")
    print("sigma_k^2 = (-2)^2 * P(k, k-2) + (-1)^2 * P(k, k-1) + (1)^2 * P(k, k+1) + (2)^2 * P(k, k+2)")
    print("sigma_k^2 = 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)")
    print("sigma_k^2 = 1 + 1/4 - c/k + 1/4 + c/k + 1")
    sigma_sq = 2.5
    print(f"After simplification, the expected squared jump is: sigma_k^2 = {sigma_sq}")
    print("-" * 50)

    # Step 3: Apply the criterion for transience
    print("Step 3: Apply the Lamperti-Pakes criterion for transience.")
    print("A chain is transient if lim_{k->inf} (k * mu_k) > (1/2) * lim_{k->inf} (sigma_k^2).")
    print("We calculate the limits:")
    print("lim_{k->inf} (k * mu_k) = lim_{k->inf} (k * 2*c/k) = 2*c")
    print(f"lim_{k->inf} (sigma_k^2) = {sigma_sq} (since it's constant for large k)")
    print("-" * 50)

    # Step 4: Solve the inequality for c
    print("Step 4: Substitute the limits into the criterion and solve for c.")
    rhs_val = 0.5 * sigma_sq
    # The final equation is 2*c > rhs_val
    print(f"The inequality for transience is: 2 * c > (1/2) * {sigma_sq}")
    print(f"Which simplifies to the final equation: 2 * c > {rhs_val}")
    
    c_inf_val = rhs_val / 2
    print(f"Solving for c gives: c > {rhs_val} / 2")
    print(f"c > {c_inf_val}")
    
    c_inf_frac = fractions.Fraction(c_inf_val).limit_denominator()
    print(f"In fractional form, the condition is: c > {c_inf_frac.numerator}/{c_inf_frac.denominator}")
    print("-" * 50)

    # Step 5: State the conclusion
    print("Step 5: Determine the infimum.")
    print(f"The chain is transient if c > {c_inf_frac}.")
    print(f"The set of values for c where the chain is transient is ({c_inf_frac}, infinity).")
    print(f"The infimum (greatest lower bound) of this set is the critical value itself.")
    
    final_answer = c_inf_frac
    print(f"\nThe infimum is {final_answer}.")

find_infimum_c()
<<<5/8>>>