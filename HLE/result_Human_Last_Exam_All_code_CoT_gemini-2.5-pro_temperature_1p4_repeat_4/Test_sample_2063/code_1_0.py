def solve_quantum_dna_probability():
    """
    Calculates the limiting probability P(n) as n -> infinity for the
    quantum DNA replication process.
    """

    # 1. Define the parameters of the problem.
    # There are 8 distinct nucleotide bases.
    num_bases = 8
    
    print(f"Step 1: The system has {num_bases} distinct bases (from 0 to {num_bases - 1}).")

    # 2. Define the failure and success conditions for a single base.
    # The process fails for a base 'x' if its count c_x satisfies: c_x mod 8 = x.
    # As n -> infinity, the probability of this specific failure condition for any
    # single base approaches 1/8, as the count becomes uniformly distributed mod 8.
    prob_failure_single_base = 1 / num_bases
    
    print(f"Step 2: The process fails for a single base 'x' if its count `c_x` satisfies `c_x mod {num_bases} = x`.")
    print(f"Step 3: As the sequence length `n` approaches infinity, the probability of this failure for any single base approaches 1/{num_bases}.")

    # 3. The probability of success for a single base is the complement.
    prob_success_single_base = 1 - prob_failure_single_base
    
    print(f"Step 4: Therefore, the probability of success for a single base (i.e., `c_x mod {num_bases} != x`) is 1 - 1/{num_bases} = {int(prob_success_single_base * num_bases)}/{num_bases}.")
    
    # 4. For total success, this condition must hold for all 8 bases.
    # In the limit n -> infinity, these events become statistically independent.
    # So, the total probability is the product of individual probabilities.
    limit_probability = prob_success_single_base ** num_bases
    
    print(f"Step 5: For the entire sequence to be replicated successfully, this condition must hold for all {num_bases} bases simultaneously.")
    print("Step 6: In the limit, these events are independent, so we multiply their probabilities.")

    # 5. Calculate the final numbers for the expression.
    numerator = (num_bases - 1)**num_bases
    denominator = num_bases**num_bases

    print("\nFinal Calculation:")
    print(f"The final probability is the success rate for one base raised to the power of the number of bases.")
    # The prompt requires showing each number in the final equation.
    final_equation_str = f"P(limit) = ({num_bases - 1}/{num_bases})^{num_bases} = {num_bases - 1}^{num_bases} / {num_bases}^{num_bases}"
    print(final_equation_str)
    
    result_str = f"P(limit) = {numerator} / {denominator}"
    print(result_str)
    
    final_value_str = f"P(limit) = {limit_probability}"
    print(final_value_str)

solve_quantum_dna_probability()