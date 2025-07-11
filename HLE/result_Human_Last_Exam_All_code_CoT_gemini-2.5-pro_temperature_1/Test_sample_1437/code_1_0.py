import sympy as sp

def solve_lindhard_function():
    """
    Calculates the universal value of the normalized Lindhard polarization function
    at zero frequency and zero momentum transfer for a 3D electron gas at T=0.
    """
    # Explain the context of the calculation
    print("In many-body physics, the static Lindhard function Pi_0(k, 0) describes the screening properties of an electron gas.")
    print("In the limit of zero momentum transfer (k -> 0), its value is Pi_0(0, 0) = -g(epsilon_F), where g(epsilon_F) is the density of states at the Fermi level.")
    print("This value depends on the material's properties. However, the *normalized* Lindhard function, Chi_0 = Pi_0 / g(epsilon_F), has a universal value in this limit.")
    print("\nWe will calculate this universal value.")
    print("-" * 50)

    # The static Lindhard function can be written as Pi_0(k, 0) = -g(epsilon_F) * F(x),
    # where x = k / (2*k_F) is a dimensionless variable.
    # The dimensionless function F(x) is:
    # F(x) = 1/2 + (1 - x**2) / (4*x) * ln| (1+x) / (1-x) |
    # We need to find the limit of F(x) as x -> 0.

    # Define the symbolic variable x
    x = sp.Symbol('x')

    # Define the dimensionless function F(x)
    # Note: sp.log is the natural logarithm (ln)
    F_x = sp.Rational(1, 2) + (1 - x**2) / (4 * x) * sp.log(sp.Abs((1 + x) / (1 - x)))

    print("Step 1: Define the dimensionless part of the Lindhard function F(x).")
    print(f"F(x) = {F_x}")

    # Calculate the limit of F(x) as x approaches 0
    limit_F_at_zero = sp.limit(F_x, x, 0)

    print("\nStep 2: Calculate the limit of F(x) as x -> 0 (which corresponds to k -> 0).")
    print(f"The limit is: {limit_F_at_zero}")

    # The normalized Lindhard function at k=0 is Chi_0(0,0) = -F(0)
    final_value = -limit_F_at_zero

    print("\nStep 3: The normalized Lindhard function Chi_0(k=0, w=0) is equal to -F(0).")
    print("\nFinal Equation:")
    print(f"Chi_0(k=0, w=0) = - ( {limit_F_at_zero} ) = {final_value}")

solve_lindhard_function()