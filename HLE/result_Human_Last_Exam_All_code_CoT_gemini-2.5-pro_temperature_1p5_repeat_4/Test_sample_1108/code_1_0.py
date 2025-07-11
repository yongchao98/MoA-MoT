import math

def solve_linearization():
    """
    Calculates the coefficients of the linearized Leslie's predator-prey system
    at its non-trivial equilibrium point.
    """
    # Step 1: Define system parameters
    a, b, c, h, m = 2, 1, 1, 1, 1

    # Step 2: Find the non-trivial equilibrium point (S_eq, F_eq)
    # From dS/dt = S(h - m*S/F) = 0, we get h = m*S/F => 1 = S/F => S = F.
    # From dF/dt = F(a - b*F - c*S) = 0, we get a = b*F + c*S => 2 = F + S.
    # Substituting S = F yields 2 = F + F => 2F = 2 => F_eq = 1.
    # Since S = F, S_eq = 1.
    S_eq = 1
    F_eq = 1

    # Step 3: Define the functions for the partial derivatives of the Jacobian matrix.
    # System: f(S, F) = S(h - mS/F) = hS - mS^2/F
    #         g(S, F) = F(a - bF - cS) = aF - bF^2 - cSF
    # Partial derivatives:
    # df/dS = h - 2*m*S/F
    # df/dF = m*S^2/F^2
    # dg/dS = -c*F
    # dg/dF = a - 2*b*F - c*S

    # Step 4: Evaluate partial derivatives at the equilibrium point (S_eq, F_eq)
    # These are the coefficients a_ij.
    a11 = h - 2 * m * S_eq / F_eq
    a12 = m * (S_eq**2) / (F_eq**2)
    a21 = -c * F_eq
    a22 = a - 2 * b * F_eq - c * S_eq

    # Step 5: For linearization at an equilibrium point, the constant terms are zero.
    b11 = 0
    b22 = 0

    # Step 6: Print the list of coefficients
    print("The coefficients of the linearized system are:")
    print(f"a_11 = {a11}")
    print(f"a_12 = {a12}")
    print(f"a_21 = {a21}")
    print(f"a_22 = {a22}")
    print(f"b_11 = {b11}")
    print(f"b_22 = {b22}")

    # As requested, also show the final equations with the numbers.
    print("\nThe final linearized system in equation form is:")
    # Using integer conversion for cleaner output as all values are integers.
    print(f"x'(t) = {int(a11)}*x(t) + {int(a12)}*y(t) + {int(b11)}")
    print(f"y'(t) = {int(a21)}*x(t) + {int(a22)}*y(t) + {int(b22)}")

if __name__ == '__main__':
    solve_linearization()