import sympy
from sympy import symbols, integrate, oo, exp, gamma, Rational

def solve_decay_rate():
    """
    Calculates the time-decay exponent for ||∂xρ(t)||_L2.

    The decay rate is determined by analyzing the long-time behavior of the
    linearized Stokes transport system. For smooth initial data, the decay
    is dominated by the low-wavenumber modes.

    The squared norm ||∂xρ(t)||_L2^2 is proportional to an integral of the form:
    I(t) = ∫[0, ∞] k^2 * exp(-C * t * k^4) dk

    This script computes the integral to find its dependence on t, and from
    that, determines the decay exponent for ||∂xρ(t)||_L2.
    """
    # Define symbolic variables for the calculation.
    # k: wave number
    # t: time
    # C: a positive constant
    k, t, C = symbols('k t C', positive=True)

    # Define the integrand based on the dominant low-wavenumber behavior.
    integrand = k**2 * exp(-C * t * k**4)

    # Compute the definite integral from 0 to infinity.
    # This gives the time dependence of the squared norm ||∂xρ(t)||_L2^2.
    integral_result = integrate(integrand, (k, 0, oo))

    # The result of the integral is gamma(3/4) / (4 * C^(3/4) * t^(3/4)).
    # We are interested in the exponent of t.
    # From the result, we can see that ||∂xρ(t)||_L2^2 is proportional to t^(-3/4).
    exponent_squared_norm = Rational(-3, 4)

    # The decay exponent for the norm ||∂xρ(t)||_L2 is half of the exponent
    # for the squared norm.
    final_exponent = exponent_squared_norm / 2
    
    numerator = final_exponent.p
    denominator = final_exponent.q

    print("Step 1: The long-time decay of ||∂xρ(t)||^2 is estimated by the integral:")
    print(f"    I(t) = ∫₀^∞ k² * exp(-C*t*k⁴) dk")
    print("\nStep 2: Using symbolic integration, the result is:")
    print(f"    I(t) = {integral_result}")
    print(f"This shows that ||∂xρ(t)||^2 is proportional to t^({exponent_squared_norm}).")
    
    print("\nStep 3: The decay for the L2-norm ||∂xρ(t)|| is the square root of the decay for its square.")
    print(f"    Exponent for ||∂xρ(t)|| = ({exponent_squared_norm}) / 2 = {final_exponent}")

    print("\nConclusion:")
    print("The best time-decay for ||∂xρ(t)||_L2 follows the relation:")
    print(f"    ||∂xρ(t)||_L2  ~  t^({numerator}/{denominator})")

if __name__ == '__main__':
    solve_decay_rate()