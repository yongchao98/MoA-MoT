def solve_fourier_restriction_problem():
    """
    Calculates the largest possible value of p based on the Fourier restriction theorem for the moment curve.
    """
    # The dimension of the space R^n
    n = 3

    # The formula for the critical value of p is n * (n + 1) / 2.
    # We will calculate this step by step.
    
    n_plus_1 = n + 1
    numerator = n * n_plus_1
    denominator = 2
    p_critical = numerator / denominator

    # The final output explains the origin of the formula and shows the calculation.
    print(f"""The problem is to find the largest value of p for which no non-zero L^p function on R^3 can have its Fourier support on the moment curve.
This is a known result from Fourier restriction theory.

The critical value of p for the moment curve in R^n is given by the formula:
p = n * (n + 1) / 2

In this specific problem, the space is R^3, so we set n = {n}.

We substitute this value into the formula to find the largest possible value of p:
p = {n} * ({n} + 1) / 2
p = {n} * {n_plus_1} / 2
p = {numerator} / {denominator}
p = {int(p_critical)}

For any p > {int(p_critical)}, it is possible to construct a non-zero L^p function whose Fourier transform is supported on the curve.
For p <= {int(p_critical)}, any L^p function with its Fourier transform supported on the curve must be identically zero.

Thus, the largest possible value of p is {int(p_critical)}.
""")

solve_fourier_restriction_problem()