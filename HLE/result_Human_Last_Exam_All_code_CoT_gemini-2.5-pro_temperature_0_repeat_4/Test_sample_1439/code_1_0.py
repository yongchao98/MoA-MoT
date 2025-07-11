def find_exponent_order():
    """
    This function explains the logic for finding the order in the coupling 'u'
    at which the critical exponent nu receives its first correction in phi^4 theory.
    """
    # In mean-field theory (equivalent to the Gaussian fixed point where coupling u=0),
    # the critical exponent nu has a classical value.
    nu_0 = 0.5
    print(f"The mean-field value of the exponent (at u=0) is nu_0 = {nu_0}")

    # The exponent nu is related to the anomalous dimension of the mass operator, gamma_r(u).
    # The exact relation is: nu = 1 / (2 - gamma_r(u))
    print("\nAway from the mean-field limit, nu is given by the formula:")
    print("nu = 1 / (2 - gamma_r(u))")

    # Perturbation theory (loop expansion) shows that the leading term of gamma_r(u)
    # is linear in the coupling constant u.
    # gamma_r(u) = -C * u + O(u^2), where C is a positive constant.
    print("\nThe one-loop calculation for the anomalous dimension gives:")
    print("gamma_r(u) = -C * u^1 + (higher order terms in u)")
    print("where C is a positive constant.")

    # Substituting this back into the formula for nu:
    # nu(u) = 1 / (2 - (-C*u)) = 1 / (2 + C*u)
    print("\nSubstituting this into the formula for nu gives:")
    print("nu(u) = 1 / (2 + C*u^1)")

    # To find the first correction, we can Taylor expand this expression for small u.
    # nu(u) = (1/2) * (1 / (1 + C/2 * u)) ≈ (1/2) * (1 - C/2 * u) = 1/2 - C/4 * u
    print("\nExpanding for small u, we get the final equation for nu:")
    # The prompt requires printing each number in the final equation.
    # The equation is nu ≈ 1/2 - (C/4) * u^1
    print("nu ≈ (1 / 2) - (C / 4) * u^1")

    # The first correction term is -(C/4)*u, which is of order u^1.
    order = 1
    print(f"\nThe first correction to the mean-field value of 1/2 is proportional to u^{order}.")
    print(f"Therefore, the order of the initial non-vanishing contribution is {order}.")

find_exponent_order()