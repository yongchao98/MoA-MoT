import numpy as np
import sympy

def solve_and_explain():
    """
    This script determines a constant lower bound for d(t,x) based on its initial condition.
    """
    print("Step 1: Deriving the evolution equation for d = u_x.")
    print("The PDE is: u_t + (u(1-u)exp(-u_bar))_x = 0, where u_bar is the integral of u.")
    print("Differentiating with respect to x gives an equation for d = u_x of the form:")
    print("d_t + ( ... )*d_x = G(u, d) * exp(-u_bar)")
    print("Along characteristics, D'(t) = exp(-u_bar) * G(u, D(t)), where D is the value of d.")
    print("After calculation, G(u, D) = 2*D**2 + (u - u**2)*D + u**3 - u**4.\n")

    print("Step 2: Applying the comparison principle.")
    print("The initial minimum of d is -0.5. We check if d can decrease below this value.")
    print("This depends on the sign of G(u, d) when d = -0.5.")
    print("Let's analyze h(u) = G(u, -0.5) for u in [0, 1].")

    u = sympy.Symbol('u')
    d_val = -0.5
    h_sympy = 2*d_val**2 + (u - u**2)*d_val + u**3 - u**4
    print(f"h(u) = G(u, {d_val}) = {h_sympy}\n")

    print("Step 3: Finding the minimum of h(u) on the interval [0, 1].")
    h_prime_sympy = sympy.diff(h_sympy, u)
    print(f"The derivative h'(u) is: {h_prime_sympy}")

    # The coefficients of h'(u) = -4*u**3 + 3*u**2 + 1*u - 0.5
    coeffs = [-4, 3, 1, -0.5]
    roots = np.roots(coeffs)
    print(f"The roots of the derivative h'(u) = 0 are approximately: {[np.round(r, 4) for r in roots]}")

    # Filter for real roots within the interval [0, 1]
    critical_points_in_interval = [r.real for r in roots if np.isreal(r) and 0 <= r.real <= 1]
    print(f"The critical points within the interval [0, 1] are: {[np.round(p, 4) for p in critical_points_in_interval]}\n")

    # Lambdify the symbolic expression for numerical evaluation
    h_func = sympy.lambdify(u, h_sympy, 'numpy')

    # Points to check for minimum value: boundaries and critical points
    points_to_check = [0.0] + critical_points_in_interval + [1.0]
    values_at_points = h_func(points_to_check)

    print("Evaluating h(u) at the critical points and at the boundaries u=0, u=1:")
    for point, value in zip(points_to_check, values_at_points):
        print(f"h({point:.4f}) = {value:.4f}")

    min_value = np.min(values_at_points)
    print(f"\nThe minimum value of h(u) on [0, 1] is approximately {min_value:.4f}.\n")

    print("Step 4: Conclusion.")
    if min_value >= 0:
        print("Since the minimum value of G(u, -0.5) is non-negative, the time derivative of d is always non-negative whenever d = -0.5.")
        print("This means that if the minimum of d(t,x) reaches -0.5, it can only increase or stay constant.")
        print("Since the initial minimum is -0.5, it cannot drop below this value.")
        lower_bound = -0.5
        print("\nA constant lower bound is therefore the initial minimum.")
        print(f"The constant lower bound is: {lower_bound}")
    else:
        print("A lower bound of -0.5 cannot be guaranteed by this method.")

if __name__ == '__main__':
    solve_and_explain()