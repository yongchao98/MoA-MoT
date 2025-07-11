import scipy.integrate
import math

def solve():
    """
    Calculates the exact value of f_Z(0.2) by solving the derived integral expression.
    """

    # G(y|x) is the CDF of the distance D = |X - x| where X ~ U[0,1].
    # It's the length of the intersection of the interval [x-y, x+y] and [0,1].
    def G(y, x):
        if y < 0:
            return 0
        return max(0, min(x + y, 1) - max(x - y, 0))

    # g(y|x) is the PDF of the distance D, which is the derivative of G(y|x) w.r.t. y.
    def g(y, x):
        if y < 0:
            return 0
        # The derivative is the sum of the derivatives of the two terms in G.
        # d/dy(min(x+y, 1)) is 1 if x+y < 1, and 0 otherwise.
        # d/dy(-max(x-y, 0)) is 1 if x-y > 0, and 0 otherwise.
        val = 0
        if x + y < 1:
            val += 1
        if x - y > 0:
            val += 1
        return val

    # f_Y_given_x is the conditional PDF of Y=D_(2) given X1=x.
    # Y is the second order statistic of 3 i.i.d. distances D.
    # The general formula for the k-th order statistic of n variables is n!/((k-1)!(n-k)!) * F^(k-1) * (1-F)^(n-k) * f
    # For n=3, k=2, this is 6 * F * (1-F) * f.
    def f_Y_given_x(y, x):
        if y <= 0:
            return 0
        G_val = G(y, x)
        g_val = g(y, x)
        return 6 * G_val * (1 - G_val) * g_val

    # This is the integrand of the inner integral for f_Z(z).
    def inner_integrand(x, y, z):
        return f_Y_given_x(y, x)

    # This is the integrand of the outer integral for f_Z(z).
    def outer_integrand(y, z):
        # The integrand is I(y)/y. We handle y=0 case.
        if y <= 1e-9: # Avoid division by zero
            return 0
        
        # Define the integration limits for the inner integral over x.
        lower_x = max(0, z - y)
        upper_x = min(1, z + y)
        
        # If the integration range is empty, the integral is zero.
        if lower_x >= upper_x:
            return 0
        
        # Perform the inner integration numerically.
        inner_integral_val, _ = scipy.integrate.quad(
            inner_integrand, lower_x, upper_x, args=(y, z)
        )
        
        return inner_integral_val / y

    # The value of z at which to evaluate the PDF.
    z_val = 0.2
    
    # Perform the outer integration to get the final result.
    # We set a higher limit for subdivisions to handle the complexity of the integrand.
    result, _ = scipy.integrate.quad(outer_integrand, 0, 1, args=(z_val,), limit=100)
    
    # The final equation is f_Z(z) = integral from 0 to 1 of (1/y * integral from max(0,z-y) to min(1,z+y) of f_Y(y|x) dx) dy
    # We are calculating this for z = 0.2
    print(f"The value of the PDF f_Z(z) at z = {z_val} is calculated based on the equation:")
    print("f_Z(0.2) = \u222B_0^1 (1/y) * (\u222B_{max(0, 0.2-y)}^{min(1, 0.2+y)} f_Y(y|x) dx) dy")
    print(f"\nThe calculated value is: {result}")

solve()