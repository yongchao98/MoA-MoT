import math

def prove_moment_problem():
    """
    This function prints a proof for the following problem:

    Suppose f: R -> R is a Schwartz class function such that the integral of x^k * f(x) over R is 0 for all non-negative integers k.
    Does it follow that f = 0?

    The function will print the logical steps to conclude that the answer is yes.
    """

    print("The statement is: If f is a Schwartz function such that ∫x^k * f(x) dx = 0 for all k ∈ {0, 1, 2, ...}, then f(x) = 0 for all x.")
    print("The answer is YES. Here is the proof:\n")

    print("--- Step 1: The integral against any polynomial is zero. ---")
    print("We are given that ∫x^k * f(x) dx = 0 for all non-negative integers k.")
    print("Any polynomial P(x) is a finite linear combination of monomials x^k:")
    print("P(x) = a_n*x^n + ... + a_1*x + a_0*x^0.")
    print("By the linearity of the integral, we can write:")
    print("∫P(x)f(x)dx = ∫(a_n*x^n + ... + a_0)f(x)dx = a_n*∫x^n*f(x)dx + ... + a_0*∫x^0*f(x)dx.")
    print("Since each term ∫x^k*f(x)dx is given to be 0, the entire sum is 0.")
    
    # Example equation with numbers as requested by the prompt format
    a2, a0, k2_moment, k0_moment, result = 7, 4, 0, 0, 0
    print(f"\nFor example, if P(x) = {a2}x^2 + {a0}, the integral is:")
    print(f"∫({a2}x^2 + {a0})f(x)dx = {a2} * ∫x^2*f(x)dx + {a0} * ∫x^0*f(x)dx = {a2}*({k2_moment}) + {a0}*({k0_moment}) = {result}.\n")
    print("So, we have established that ∫P(x)f(x)dx = 0 for any polynomial P(x).\n")

    print("--- Step 2: Extending from polynomials to continuous functions with compact support. ---")
    print("Let g(x) be any continuous function with compact support. Let's say the support is contained in [-A, A].")
    print("The Weierstrass Approximation Theorem states that g(x) can be uniformly approximated by a polynomial P(x) on [-A, A].")
    print("This means for any ε > 0, there exists a polynomial P(x) such that |g(x) - P(x)| < ε for all x in [-A, A].")
    print("\nNow let's consider the integral ∫g(x)f(x)dx:")
    print("∫g(x)f(x)dx = ∫(g(x) - P(x) + P(x))f(x)dx = ∫(g(x) - P(x))f(x)dx + ∫P(x)f(x)dx.")
    print("From Step 1, we know the second term ∫P(x)f(x)dx is 0.")
    print("So, |∫g(x)f(x)dx| = |∫(g(x) - P(x))f(x)dx|.")
    print("\nBecause g(x) has support in [-A, A], this integral is over [-A, A].")
    print("|∫g(x)f(x)dx| ≤ ∫[-A,A] |g(x) - P(x)|*|f(x)|dx < ∫[-A,A] ε*|f(x)|dx = ε * ∫[-A,A] |f(x)|dx.")
    print("Since f is a Schwartz function, it is integrable, so ∫|f(x)|dx is a finite constant.")
    print("Because we can make ε arbitrarily small, the only non-negative value that |∫g(x)f(x)dx| can be is 0.")
    print("Thus, ∫g(x)f(x)dx = 0 for any continuous function g(x) with compact support.\n")

    print("--- Step 3: Proving f(x) must be the zero function. ---")
    print("We will use proof by contradiction. Assume f(x) is not identically zero.")
    print("This means there must be some point x_0 where f(x_0) ≠ 0. Let's assume f(x_0) = c > 0.")
    print("Since f is a Schwartz function, it is continuous. By continuity, if f(x_0) > 0, there must be a small interval I = (x_0 - δ, x_0 + δ) around x_0 where f(x) > 0 for all x in I.")
    print("\nNow, we can construct a continuous function g(x) with compact support, which is entirely contained in the interval I, and where g(x) ≥ 0, but g(x) is not identically zero (e.g., a 'bump' or 'tent' function).")
    print("Consider the integral of the product g(x)f(x).")
    print("Inside the interval I, f(x) > 0 and g(x) ≥ 0 (and g(x) > 0 on some subinterval).")
    print("Outside the interval I, g(x) = 0.")
    print("Therefore, the product g(x)f(x) is non-negative everywhere, and strictly positive on some subinterval.")
    print("The integral of a non-negative, continuous function which is not identically zero must be strictly positive. So, ∫g(x)f(x)dx > 0.")
    print("\nThis is a contradiction! In Step 2, we proved that ∫g(x)f(x)dx = 0 for ALL continuous functions g(x) with compact support. Our constructed bump function is one such g(x).")
    print("The contradiction arises from our assumption that f(x) is not identically zero.")
    print("Therefore, the assumption must be false.\n")

    print("--- Conclusion ---")
    print("The only possibility is that f(x) = 0 for all x.")

if __name__ == '__main__':
    prove_moment_problem()
