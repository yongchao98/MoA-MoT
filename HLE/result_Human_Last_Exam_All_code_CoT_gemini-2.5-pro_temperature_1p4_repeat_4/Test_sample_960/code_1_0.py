from fractions import Fraction

def calculate_probability_bound():
    """
    This function calculates the upper bound for the probability that the production
    process reaches a point where exactly 50% of the products are good.

    The probability is calculated by summing the probabilities of stopping at each
    possible time T. It can be shown that T must be odd (T=2k-1) and the
    probability of stopping at T=2k-1 is P(T=2k-1) = 1 / ((2k-1)*(2k+1)).

    The total probability is the sum of this series from k=1 to infinity, which
    is a telescoping series that converges to 1/2.
    """
    print("Calculating the probability by summing the series P(T < infinity) = sum_{k=1 to inf} P(T=2k-1).")
    print("The general term for the probability of stopping at time T=2k-1 is 1/((2k-1)*(2k+1)).\n")

    num_terms = 1000  # Number of terms for the numerical approximation
    total_prob = Fraction(0)
    
    equation_terms = []
    
    for k in range(1, num_terms + 1):
        term = Fraction(1, (2 * k - 1) * (2 * k + 1))
        total_prob += term
        if k <= 5: # Store first 5 terms for the equation string
            equation_terms.append(f"{term.numerator}/{term.denominator}")

    print(f"The sum of the first {num_terms} terms is: {total_prob.numerator}/{total_prob.denominator}")
    print(f"This is approximately: {float(total_prob)}")
    print("The exact sum of the infinite series is 1/2, which is 0.5.\n")
    print("The upper bound for the probability is its exact value, 1/2.")
    
    # Print the final equation with the first few terms as requested
    final_equation_str = " + ".join(equation_terms)
    final_equation_str += " + ... = 1/2"
    
    print("\nThe final equation is:")
    print(final_equation_str)


calculate_probability_bound()