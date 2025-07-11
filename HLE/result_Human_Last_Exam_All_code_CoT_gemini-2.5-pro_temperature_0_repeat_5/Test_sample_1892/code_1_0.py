import math

def solve_asymptotic_parameters():
    """
    Calculates the parameters alpha and beta for the asymptotic formula.

    The asymptotic formula is |A(X)| ~ c * X^alpha * (log X)^beta.
    The theory of Dirichlet series tells us that alpha = 1 and beta = k - 1,
    where k is the order of the pole at s=1 of the Dirichlet series
    D(s) = sum(f(q) / q^s), with f(q) being the number of primitive
    characters mod q whose order divides 12.

    The value of k is the average of the function g(a) = gcd(12, a-1) - 1
    over the group of units modulo 12.
    """
    m = 12

    # Find phi(m), the number of integers a such that 1 <= a < m and gcd(a, m) = 1.
    phi_m = 0
    units = []
    for a in range(1, m):
        if math.gcd(a, m) == 1:
            phi_m += 1
            units.append(a)

    # Calculate the sum of g(a) for a in the group of units.
    # Note: gcd(m, 0) is defined as m in this context.
    term_sum = 0
    print("Calculating the order of the pole, k:")
    print(f"m = {m}, phi(m) = {phi_m}")
    print(f"Units mod {m}: {units}")
    print("-" * 20)
    for a in units:
        # The term is gcd(12, a-1) - 1.
        # For a=1, a-1=0. math.gcd(12, 0) is 12.
        term = math.gcd(m, a - 1) - 1
        term_sum += term
        print(f"For a = {a:2d}, term = gcd({m}, {a-1:2d}) - 1 = {term+1:2d} - 1 = {term:2d}")

    # The order of the pole k is the average of these terms.
    k = term_sum / phi_m
    print("-" * 20)
    print(f"Sum of terms = {term_sum}")
    print(f"k = Sum / phi(m) = {term_sum} / {phi_m} = {k}")
    print("-" * 20)

    # Determine alpha and beta
    alpha = 1
    beta = k - 1

    print("The asymptotic formula is of the form c * X^alpha * (log X)^beta")
    print(f"alpha = {int(alpha)}")
    print(f"beta = k - 1 = {int(k)} - 1 = {int(beta)}")
    print("-" * 20)

    # The final result is the sum of alpha and beta.
    result = alpha + beta
    print(f"The sum alpha + beta is {int(alpha)} + {int(beta)} = {int(result)}")

solve_asymptotic_parameters()