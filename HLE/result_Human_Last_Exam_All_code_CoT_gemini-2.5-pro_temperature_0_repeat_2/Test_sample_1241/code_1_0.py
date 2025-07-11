def solve_steady_state_probability():
    """
    Calculates the sum of steady-state probabilities P0(+inf) + P1(+inf)
    for the given Kolmogorov-Chepmen system.
    """
    # Given parameters
    lambda_01 = 0.019
    lambda_10 = 0.65
    
    # The problem simplifies to solving a system of linear equations for the
    # steady-state probabilities p_i = P_i(+inf).
    # The derivation leads to the following expression for p0 + p1:
    # p0 + p1 = ((lambda_10 / lambda_01) + 1) / ((lambda_10 / lambda_01) + 3)

    # Calculate the ratio lambda_10 / lambda_01
    ratio = lambda_10 / lambda_01

    # Calculate the final result for p0 + p1
    result = (ratio + 1) / (ratio + 3)

    # Print the final equation with the numbers and the result
    print(f"The final equation for P0(+inf) + P1(+inf) is derived from the steady-state conditions.")
    print(f"The calculation is as follows:")
    print(f"P0(+inf) + P1(+inf) = (({lambda_10} / {lambda_01}) + 1) / (({lambda_10} / {lambda_01}) + 3) = {result}")

solve_steady_state_probability()