import numpy as np

def analyze_lower_bound():
    """
    This function analyzes the polynomial P(u) to verify that -1 is a lower bound.
    The condition for C=-1 to be a lower bound is that the polynomial
    P(u) = u^4 - u^3 - 5u^2 + 3u + 2 is non-negative for u in [0, 1].
    """

    # Define the polynomial P(u) and its coefficients
    # P(u) = 1*u^4 - 1*u^3 - 5*u^2 + 3*u + 2
    p_coeffs = [1, -1, -5, 3, 2]
    P = np.poly1d(p_coeffs)

    # We argue that P(u) is concave on [0,1], so its minimum is at the endpoints.
    # To show this, we analyze the second derivative P''(u).
    P_double_prime = P.deriv(2)
    # P''(u) = 12u^2 - 6u - 10

    print("--- Analysis of the lower bound C = -1 ---")
    print(f"We need to check if the polynomial P(u) = {P}")
    print("is non-negative for u in the interval [0, 1].")
    print("\nStep 1: Analyze the concavity of P(u) by checking the sign of its second derivative:")
    print(f"P''(u) = {P_double_prime}")
    
    # Check values of P''(u) at endpoints of [0, 1]
    p_double_prime_at_0 = P_double_prime(0)
    p_double_prime_at_1 = P_double_prime(1)
    
    # Roots of P''(u) = 12u^2 - 6u - 10 = 0 are (6 +/- sqrt(36 + 480))/24
    # Positive root is (6 + sqrt(516))/24 > 1. Negative root is < 0.
    # So P''(u) does not change sign on [0, 1].
    print(f"Value of P''(0) = {p_double_prime_at_0}")
    print(f"Value of P''(1) = {p_double_prime_at_1}")
    print("Since P''(u) is a parabola opening upwards and its roots are outside [0, 1],")
    print("and P''(0) and P''(1) are negative, P''(u) is negative for all u in [0, 1].")
    print("This means P(u) is a concave function on [0, 1].")
    
    print("\nStep 2: Find the minimum of the concave function P(u) on [0, 1].")
    print("The minimum of a concave function on a closed interval occurs at its endpoints.")
    
    p_at_0 = P(0)
    p_at_1 = P(1)
    
    print(f"Value of P(u) at u = 0: P(0) = {int(p_at_0)}")
    print(f"Value of P(u) at u = 1: 1^4 - 1^3 - 5(1)^2 + 3(1) + 2 = {int(p_at_1)}")

    min_p_on_interval = min(p_at_0, p_at_1)
    
    print(f"\nThe minimum value of P(u) on [0, 1] is {min_p_on_interval}.")
    
    if min_p_on_interval >= 0:
        print("Since the minimum value is non-negative, P(u) >= 0 for all u in [0, 1].")
        print("This confirms that C = -1 is a valid lower bound for d(t,x).")
        
    final_answer = -1
    print(f"\nThe constant lower bound for d(t,x) is {final_answer}.")

if __name__ == '__main__':
    analyze_lower_bound()
