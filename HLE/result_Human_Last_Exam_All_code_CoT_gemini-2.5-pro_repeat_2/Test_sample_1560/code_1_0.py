import numpy as np

def solve_ququint_problem():
    """
    Calculates the final state and measurement probabilities for the given ququint system.
    """
    # The coefficients of the final state are derived by applying the gate Q to the initial state |ψ>.
    # Q|ψ> = Q * (1/√5) * (|0> + |1> + |2> + |3> + |4>)
    #      = (1/√5) * (Q|0> + Q|1> + Q|2> + Q|3> + Q|4>)
    #      = (1/√10) * [(|1>+|2>) + (|0>+|3>) + (|1>+|4>) + (|2>+|0>) + (|3>+|2>)]
    #      = (1/√10) * [2|0> + 2|1> + 3|2> + 2|3> + 1|4>]
    # We can work with the unnormalized coefficients [2, 2, 3, 2, 1] and normalize at the end.
    
    coeffs = np.array([2, 2, 3, 2, 1])
    
    # The square of the norm is the sum of the squares of the coefficients.
    # This is used for the normalization factor 1/sqrt(norm_sq).
    norm_sq = np.sum(coeffs**2)
    
    print("The final state of the system after applying gate Q is |ψ'>.")
    print("The equation for this state, after normalization, is:")
    
    # Build and print the final state equation term by term
    equation_str = "|ψ'> = "
    terms = []
    for i, coeff in enumerate(coeffs):
        # We output each number in the equation, showing the coefficients and the normalization factor.
        terms.append(f"({coeff}/√{norm_sq})|{i}>")
    equation_str += " + ".join(terms)
    print(equation_str)
    
    print("\nThe probability of measuring the system in each basis state |i> is P(|i>) = |coefficient of |i>|².")
    
    # Calculate and print the measurement probabilities
    probabilities = (coeffs**2) / norm_sq
    for i, prob in enumerate(probabilities):
        # We output each number used in the probability calculation.
        print(f"P(|{i}>) = |{coeffs[i]}/√{norm_sq}|² = {coeffs[i]**2}/{norm_sq} ≈ {prob:.4f}")

    # Verify that the probabilities sum to 1
    print(f"\nTotal Probability = {np.sum(probabilities):.2f}")

solve_ququint_problem()