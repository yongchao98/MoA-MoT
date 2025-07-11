import sympy

def calculate_variance_of_coxeter_length():
    """
    Calculates the variance of the Coxeter length on the hyperoctahedral group B_3
    using the Poincaré polynomial and symbolic differentiation.
    """
    # 1. Define the symbol q and the rank n
    q = sympy.Symbol('q')
    n = 3

    print(f"Calculating the variance for the Coxeter length on B_{n}.")
    print("The method uses the Poincaré polynomial P(q) for B_n.")

    # 2. Define the Poincaré polynomial P(q) for B_3
    # P(q) = product_{i=1 to n} [2i]_q, where [k]_q = (1 - q^k) / (1 - q)
    poly_factors = []
    for i in range(1, n + 1):
        poly_factors.append((1 - q**(2 * i)) / (1 - q))

    P = 1
    for factor in poly_factors:
        P *= factor

    print(f"For n={n}, the polynomial is P(q) = [2]q * [4]q * [6]q.")

    # 3. Calculate the first and second derivatives of P(q)
    P_prime = sympy.diff(P, q)
    P_double_prime = sympy.diff(P_prime, q)

    # 4. Evaluate P, P', and P'' at q=1 using limits to avoid division by zero
    P_1 = sympy.limit(P, q, 1)
    P_prime_1 = sympy.limit(P_prime, q, 1)
    P_double_prime_1 = sympy.limit(P_double_prime, q, 1)

    print("\n--- Intermediate Values ---")
    print(f"The size of the group |B_{n}| is P(1) = {P_1}")
    print(f"The value of the first derivative P'(1) is: {P_prime_1}")
    print(f"The value of the second derivative P''(1) is: {P_double_prime_1}")

    # 5. Calculate mean and variance using the standard formulas
    mean_L = P_prime_1 / P_1
    
    # E[L^2] = (P''(1) + P'(1)) / P(1)
    mean_L_squared = (P_double_prime_1 + P_prime_1) / P_1
    
    # Var(L) = E[L^2] - (E[L])^2
    variance_L = mean_L_squared - mean_L**2
    
    # Convert to rational for clean printing
    mean_L_rational = sympy.Rational(mean_L)
    mean_L_squared_rational = sympy.Rational(mean_L_squared)
    variance_L_rational = sympy.Rational(variance_L)


    print("\n--- Variance Calculation ---")
    print("The formula for variance is Var(L) = E[L^2] - (E[L])^2")
    print(f"The mean E[L] = P'(1)/P(1) = {P_prime_1}/{P_1} = {mean_L_rational}")
    print(f"The mean of squares E[L^2] = (P''(1) + P'(1))/P(1) = ({P_double_prime_1} + {P_prime_1})/{P_1} = {mean_L_squared_rational}")

    # Final equation printout
    print("\nFinal Equation:")
    print(f"{mean_L_squared_rational} - ({mean_L_rational})^2 = {mean_L_squared_rational} - {mean_L_rational**2} = {variance_L_rational}")
    print(f"The variance as a floating-point number is: {float(variance_L)}")

calculate_variance_of_coxeter_length()