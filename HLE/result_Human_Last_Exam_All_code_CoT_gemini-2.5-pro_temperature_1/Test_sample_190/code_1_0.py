import fractions

def solve_markov_chain_transience():
    """
    Calculates the infimum of c for a Markov chain to be transient.

    The method involves:
    1. Calculating the mean drift mu_k, which for large k is of the form a/k.
    2. Calculating the limiting variance of the jumps, sigma^2.
    3. Applying the transience criterion: a > sigma^2 / 2.
    4. Solving the resulting inequality for c.
    """

    print("--- Step 1: Calculate the mean drift coefficient 'a' ---")
    print("The mean drift mu_k is E[X_{n+1} - k | X_n = k].")
    print("The jumps from state k are d = {-2, -1, 1, 2} with probabilities:")
    print("P(d=-2) = 1/4")
    print("P(d=-1) = 1/4 - c/k")
    print("P(d=+1) = 1/4 + c/k")
    print("P(d=+2) = 1/4")
    print("\nmu_k = sum(d * P(d))")
    print("mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print("mu_k = -1/2 - 1/4 + c/k + 1/4 + c/k + 1/2")
    print("mu_k = (2*c)/k\n")

    # The drift is of the form a/k. The coefficient 'a' is 2*c.
    a_coeff_c = 2
    print(f"From this, we see the drift is of the form a/k, where a = {a_coeff_c}*c.\n")


    print("--- Step 2: Calculate the limiting variance sigma^2 ---")
    print("sigma^2 is the limit as k->infinity of the expected squared jump, E[(X_{n+1} - k)^2 | X_n = k].")
    print("E_k[(X_{n+1}-k)^2] = sum(d^2 * P(d))")
    print("= (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)")
    print("= 4*(1/4) + 1*(1/4 - c/k) + 1*(1/4 + c/k) + 4*(1/4)")
    print("= 1 + 1/4 - c/k + 1/4 + c/k + 1")
    print("= 2 + 1/2 = 5/2\n")

    sigma_sq = fractions.Fraction(5, 2)
    print(f"The limit sigma^2 is {sigma_sq}.\n")


    print("--- Step 3: Apply the transience criterion and solve for c ---")
    print("A 1D random walk with drift a/k and variance sigma^2 is transient if: a > sigma^2 / 2")
    
    # Values for the equation
    sigma_sq_div_2 = sigma_sq / 2
    
    print("\nSubstituting our derived values 'a' and 'sigma^2':")
    # Outputting each number in the final equation
    print(f"The inequality is: {a_coeff_c}*c > {sigma_sq.numerator}/{sigma_sq.denominator} / 2")
    print(f"Which simplifies to: {a_coeff_c}*c > {sigma_sq_div_2.numerator}/{sigma_sq_div_2.denominator}")
    
    # Solve for c
    inf_c = sigma_sq_div_2 / a_coeff_c
    
    print("\nSolving for c:")
    print(f"c > ({sigma_sq_div_2.numerator}/{sigma_sq_div_2.denominator}) / {a_coeff_c}")
    print(f"c > {inf_c.numerator}/{inf_c.denominator}\n")

    print("The set of values for c where the Markov chain is transient is (5/8, infinity).")
    print("Therefore, the infimum of this set is 5/8.")
    print("\nFinal Answer:")
    print(f"{inf_c.numerator}/{inf_c.denominator}")

if __name__ == "__main__":
    solve_markov_chain_transience()
    final_answer = fractions.Fraction(5, 8)
    print(f"\n<<<{float(final_answer)}>>>")
