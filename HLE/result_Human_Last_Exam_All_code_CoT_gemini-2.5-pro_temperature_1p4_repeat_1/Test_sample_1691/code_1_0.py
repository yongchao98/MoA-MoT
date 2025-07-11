import math

def solve():
    """
    Calculates the coefficients and exponents for the asymptotic formula of the integral.
    """
    # Parameters from the integral denominator: epsilon + p0*x^n + p1*x^m + ...
    p0 = 9.0
    n = 5.0
    p1 = 5.0
    m = 6.0

    # Exponent of the leading term
    a = (1.0 - n) / n

    # Coefficient of the leading term
    # A = p0**(-1/n) * integral(dy / (1+y**n)) from 0 to inf
    # The integral is known to be (pi/n) / sin(pi/n)
    A = math.pow(p0, -1.0 / n) * (math.pi / n) / math.sin(math.pi / n)

    # Exponent of the first correction term
    # This is derived from the expansion, which gives a power of epsilon
    # as -(n-1)/n + (m-n)/n = (m - 2n + 1)/n
    b = (m - 2.0 * n + 1.0) / n

    # Coefficient of the first correction term, B.
    # The derivation shows that the coefficient is:
    # B = -p1 * p0**(-(m+1)/n) * integral(y**m / (1+y**n)**2)
    # The integral part evaluates to (1/n) * Beta((m+1)/n, 2 - (m+1)/n)
    # Using the reflection formula for the Beta function, this becomes:
    # Beta((m+1)/n, 2-(m+1)/n) = Beta(7/5, 3/5) which is (2*pi)/(5*sin(2*pi/5))
    # The term p1=5 gets canceled by the 1/n=1/5 from the integral evaluation.
    # So B = -p0**(-(m+1)/n) * Beta((m+1)/n, 2-(m+1)/n)
    
    # We use m=6, n=5. Beta arg1=(6+1)/5=1.4, arg2=2-1.4=0.6
    # We use Beta(x,y) = Gamma(x)Gamma(y)/Gamma(x+y). Here, x+y=2, Gamma(2)=1.
    # Gamma(1.4)=0.4*Gamma(0.4), Gamma(0.6)=Gamma(1-0.4).
    # Gamma(0.4)Gamma(0.6) = pi/sin(0.4*pi)
    # Beta(1.4, 0.6) = 0.4 * pi / sin(0.4*pi) = (2/5)*pi / sin(2*pi/5)
    beta_val = (2.0 * math.pi) / (n * math.sin(2.0 * math.pi / n))
    B = -math.pow(p0, -(m + 1.0) / n) * beta_val
    
    # Print the resulting analytical formula with numerical values.
    # Note that a = -4/5 = -0.8 and b = -3/5 = -0.6
    print("An approximate analytical formula for I(epsilon) for small epsilon is:")
    print(f"I(epsilon) ≈ {A:.4f} * epsilon^({a:.2f}) + {B:.4f} * epsilon^({b:.2f})")
    print("\nWhich is:")
    print(f"I(epsilon) ≈ {A:.4f} / epsilon^{abs(a):.2f} {B:.4f} / epsilon^{abs(b):.2f}")
    
    # We are asked to output the final answer in a special format.
    # Let's present the formula as the answer.
    final_answer = f"{A:.4f} * epsilon^({a:.2f}) + {B:.4f} * epsilon^({b:.2f})"
    return final_answer

final_formula = solve()
# The final answer format is just the content, not a print statement.
# To conform to the format, let's build the final string.
# Example: <<<9.8>>>
# So we output the formula string inside the angle brackets.
# This seems to be the most reasonable interpretation of the requested format.
# Let's format it nicely.
# <<<0.6885 * epsilon**(-0.80) - 0.0708 * epsilon**(-0.60)>>>
A_val = 0.6885
a_val = -0.80
B_val = -0.0708
b_val = -0.60
# This is derived in the thought process. Let's recalculate and format from the function.
import math
p0 = 9.0; n = 5.0; p1 = 5.0; m = 6.0
a = (1.0 - n) / n
A = math.pow(p0, -1.0 / n) * (math.pi / n) / math.sin(math.pi / n)
b = (m - 2.0 * n + 1.0) / n
beta_val = (2.0 * math.pi) / (n * math.sin(2.0 * math.pi / n))
B = -math.pow(p0, -(m + 1.0) / n) * beta_val

answer_string = f"{A:.4f}*epsilon^({a:.2f}) + {B:.4f}*epsilon^({b:.2f})"
# Since B is negative, it will print with a minus sign. Let's make it cleaner.
answer_string_clean = f"{A:.4f}*epsilon^({a:.2f}) - {abs(B):.4f}*epsilon^({b:.2f})"

print(f"<<<{answer_string_clean}>>>")
