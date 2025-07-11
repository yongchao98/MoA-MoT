import sympy

def analyze_lower_bound():
    """
    This script analyzes the polynomial that determines the lower bound of d(t,x).
    """
    u = sympy.Symbol('u')
    
    # Define the polynomial P(u) which needs to be non-negative for u in [0, 1].
    # This polynomial arises from evaluating the time derivative of the minimum of d(t,x)
    # at the candidate lower bound C = -1.
    coeffs = [1, -1, -5, 3, 2]
    P = coeffs[0]*u**4 + coeffs[1]*u**3 + coeffs[2]*u**2 + coeffs[3]*u + coeffs[4]

    print("To find the lower bound, we need to analyze the sign of the following polynomial for u in [0, 1]:")
    print(f"P(u) = {coeffs[0]}*u**4 + ({coeffs[1]})*u**3 + ({coeffs[2]})*u**2 + {coeffs[3]}*u + {coeffs[4]}")
    
    print("\nStep 1: Factor the polynomial.")
    P_factored = sympy.factor(P)
    print(f"The factored form of P(u) is: {P_factored}")

    # Extract factors
    # In Sympy, the factored form might be a product of terms.
    # For P(u) = (u - 1)*(u**3 - 5*u - 2), the factors are (u-1) and (u**3 - 5*u - 2)
    
    print("\nStep 2: Analyze the sign of each factor in the interval [0, 1].")
    
    factor1_str = "u - 1"
    print(f"\nAnalyzing factor: {factor1_str}")
    print("For u in [0, 1), u - 1 is negative.")
    print("For u = 1, u - 1 is zero.")
    print("Thus, for u in [0, 1], the factor (u - 1) is non-positive (<= 0).")

    factor2_str = "u**3 - 5*u - 2"
    Q = u**3 - 5*u - 2
    print(f"\nAnalyzing factor: {factor2_str}")
    
    # Analyze the derivative of the second factor Q(u)
    Q_prime = sympy.diff(Q, u)
    print(f"The derivative is Q'(u) = {Q_prime}.")
    print("For u in [0, 1], u**2 is in [0, 1], so 3*u**2 is in [0, 3].")
    print("Therefore, 3*u**2 - 5 is always negative for u in [0, 1].")
    print("Since its derivative is negative, Q(u) is a strictly decreasing function on [0, 1].")
    
    # Evaluate Q(u) at the endpoints of the interval
    Q_at_0 = Q.subs(u, 0)
    Q_at_1 = Q.subs(u, 1)
    print(f"The maximum value of Q(u) on [0, 1] occurs at u = 0, which is {Q_at_0}.")
    print(f"The minimum value of Q(u) on [0, 1] occurs at u = 1, which is {Q_at_1}.")
    print("Since Q(u) is decreasing and Q(0) is negative, Q(u) is negative for all u in [0, 1].")

    print("\nStep 3: Conclude the sign of P(u).")
    print("P(u) is the product of (u - 1) and (u**3 - 5*u - 2).")
    print("On [0, 1], (u - 1) is non-positive and (u**3 - 5*u - 2) is negative.")
    print("The product of a non-positive number and a negative number is a non-negative number.")
    print("Therefore, P(u) >= 0 for all u in [0, 1].")

    print("\nFinal Conclusion:")
    lower_bound = -1
    print("The analysis shows that when d(t,x) reaches -1, its time derivative becomes non-negative, preventing it from decreasing further.")
    print(f"Given that the initial minimum of d(t,x) is -0.5, which is greater than {lower_bound}, the minimum of d(t,x) will remain greater than or equal to {lower_bound}.")
    print(f"A constant lower bound of d(t,x) is {lower_bound}.")

if __name__ == '__main__':
    analyze_lower_bound()