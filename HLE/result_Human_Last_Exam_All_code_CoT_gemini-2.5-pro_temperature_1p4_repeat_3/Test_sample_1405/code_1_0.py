import numpy as np

def solve_traffic_flow_bound():
    """
    This function numerically verifies the analytical derivation for the lower bound
    of the vehicle density derivative in the given traffic flow model.
    
    The lower bound is determined by the minimum of the smaller root of a quadratic equation
    P(y, u) = a*y^2 + b(u)*y + c(u) = 0, where y is the derivative d(t,x) and u is the density u(t,x).
    """

    # Coefficients of the quadratic P(y, u) = a*y^2 + b*y + c = 0
    a = 2.0

    def b_coeff(u):
        """Calculates the coefficient b as a function of u."""
        return -(3 * u - 5 * u**2)

    def c_coeff(u):
        """Calculates the coefficient c as a function of u."""
        return -(u**3 * (1 - u))

    def lower_root_y1(u):
        """Calculates the lower root y1 of the quadratic equation."""
        if u == 0:
            return 0.0
        # The discriminant is b^2 - 4ac
        discriminant = b_coeff(u)**2 - 4 * a * c_coeff(u)
        # The formula for the smaller root is (-b - sqrt(discriminant)) / (2a)
        root = (-b_coeff(u) - np.sqrt(discriminant)) / (2 * a)
        return root

    # We search for the minimum of y1(u) over the valid range of u, which is [0, 1].
    u_values = np.linspace(0, 1, 2001)
    y1_values = np.array([lower_root_y1(u) for u in u_values])

    # Find the minimum value and the value of u where it occurs.
    min_y1 = np.min(y1_values)
    min_u_index = np.argmin(y1_values)
    u_at_min = u_values[min_u_index]
    
    print("This script verifies the analytical derivation of the lower bound.")
    print(f"The minimum of the lower root y1(u) is found numerically at u = {u_at_min:.3f}")
    
    # As derived analytically, the minimum occurs exactly at u=1.
    # Let's print the equation for this case.
    u_final = 1.0
    a_final = a
    b_final = b_coeff(u_final)
    c_final = c_coeff(u_final)
    
    print("\nBased on the analytical result, the minimum occurs at u = 1.")
    print("The governing quadratic equation is of the form: a*y^2 + b*y + c = 0")
    print("For u = 1.0, the coefficients are:")
    print(f"a = {a_final}")
    print(f"b = {b_final}")
    print(f"c = {c_final}")
    print("The final equation for the bound is:")
    print(f"{a_final}*y^2 + {b_final}*y + {c_final} = 0")
    
    final_bound = lower_root_y1(u_final)
    print("\nThe constant lower bound of d(t,x) is the smaller root of this equation.")
    print(f"The calculated lower bound is: {final_bound}")

if __name__ == "__main__":
    solve_traffic_flow_bound()
