import sympy
from sympy import symbols, Function, Eq, cos, sin, series, S, sqrt, cancel

def solve_rayleigh_plesset_correction():
    """
    This function uses the Poincaré-Lindstedt method to find the nonlinear
    frequency correction for the Rayleigh-Plesset equation using sympy.
    """
    # Define the small parameter, time variable, and polytropic index
    eps, tau, gamma = symbols('epsilon tau gamma', real=True, positive=True)

    # Define R_i and lambda_i as symbols
    # We will determine these order by order
    lmbda0, lmbda1, lmbda2 = symbols('lambda_0 lambda_1 lambda_2')

    # Define R_i as functions of tau
    R1 = Function('R1')(tau)
    R2 = Function('R2')(tau)
    R3 = Function('R3')(tau)

    # Define the series expansions for R(tau) and lambda = omega^2
    R_series = 1 + eps * R1 + eps**2 * R2 + eps**3 * R3
    lmbda_series = lmbda0 + eps * lmbda1 + eps**2 * lmbda2

    # Derivatives with respect to tau
    R_p = R_series.diff(tau)
    R_pp = R_series.diff(tau, 2)

    # The dimensionless Rayleigh-Plesset equation is:
    # R * R_ddot + (3/2)*R_dot^2 = R**(-3*gamma) - 1
    # where derivatives are wrt non-dimensional time t.
    # With tau = omega * t, d/dt = omega * d/dtau.
    # So, R * (omega^2 * R'') + (3/2)*(omega*R')**2 = R**(-3*gamma) - 1
    # Let lambda = omega^2.
    # lambda * (R * R'' + (3/2)*R'**2) = R**(-3*gamma) - 1
    
    # Define the LHS and RHS of the equation
    LHS = lmbda_series * (R_series * R_pp + S(3)/2 * R_p**2)
    RHS = R_series**(-3 * gamma) - 1

    # Form the full equation and expand it as a series in epsilon
    full_eq = LHS - RHS
    eq_in_eps = series(full_eq, eps, n=4).removeO()

    # --- Order epsilon^1 ---
    # Equation at O(eps)
    eq_eps1 = eq_in_eps.coeff(eps, 1)
    # The equation is lmbda0*R1'' - (-3*gamma*R1) = 0 => lmbda0*R1'' + 3*gamma*R1 = 0
    # For a non-trivial periodic solution with frequency 1 in tau, R1'' + R1 = 0.
    # This requires lmbda0 = 3*gamma.
    lmbda0_sol = 3 * gamma
    
    # From initial conditions R(0)=1+eps, R'(0)=0, we get R1(0)=1, R1'(0)=0.
    # The solution to R1'' + R1 = 0 with these ICs is R1 = cos(tau).
    R1_sol = cos(tau)

    # --- Order epsilon^2 ---
    # Equation at O(eps^2)
    eq_eps2 = eq_in_eps.coeff(eps, 2)
    # Substitute known values
    eq_eps2 = eq_eps2.subs([(lmbda0, lmbda0_sol), (R1, R1_sol)])
    # The equation is of the form: lmbda0*(R2''+R2) = Forcing
    # We collect the coefficient of cos(tau) in the forcing term (secular term) and set it to zero.
    
    # Isolate terms for the R2 ODE
    forcing_R2 = - (eq_eps2 - lmbda0_sol * R2.diff(tau, 2) - 3*gamma * R2)
    forcing_R2_trig = sympy.trigsimp(forcing_R2)
    
    # The secular term is the one with cos(tau). Its coefficient must be zero.
    secular_coeff_l1 = forcing_R2_trig.coeff(cos(tau))
    # This gives lmbda1 = 0.
    lmbda1_sol = 0

    # Now we can solve for R2.
    forcing_R2_no_l1 = forcing_R2.subs(lmbda1, lmbda1_sol)
    R2_ode_rhs = cancel(forcing_R2_no_l1 / lmbda0_sol)
    
    # Solve R2'' + R2 = R2_ode_rhs using method of undetermined coefficients
    # R2_p = C0 + C2*cos(2*tau)
    C0 = R2_ode_rhs.subs(cos(2*tau), 0)
    C2 = R2_ode_rhs.coeff(cos(2*tau)) / (1 - 2**2)
    # The homogeneous solution part comes from IC: R2(0)=0, R2'(0)=0
    # R2_h = A*cos(tau) + B*sin(tau)
    # R2(0) = A + C0 + C2 = 0 => A = -C0 - C2
    A = -C0 - C2
    B = 0 # from R2'(0)=0
    R2_sol = A * cos(tau) + B * sin(tau) + C0 + C2 * cos(2 * tau)
    
    # --- Order epsilon^3 ---
    # Equation at O(eps^3)
    eq_eps3 = eq_in_eps.coeff(eps, 3)
    
    # Substitute all known solutions
    substitutions = [
        (lmbda0, lmbda0_sol), (lmbda1, lmbda1_sol),
        (R1.diff(tau, 2), R1_sol.diff(tau, 2)), (R1.diff(tau), R1_sol.diff(tau)), (R1, R1_sol),
        (R2.diff(tau, 2), R2_sol.diff(tau, 2)), (R2.diff(tau), R2_sol.diff(tau)), (R2, R2_sol)
    ]
    eq_eps3_subbed = eq_eps3.subs(substitutions)
    
    # Isolate the forcing term for R3 ODE: lmbda0*(R3''+R3) = Forcing
    forcing_R3 = - (eq_eps3_subbed - lmbda0_sol * R3.diff(tau, 2) - 3*gamma * R3)
    
    # Expand trigonometric functions to identify the secular term (coefficient of cos(tau))
    forcing_R3_expanded = sympy.expand_trig(forcing_R3)
    
    # The secular condition is that the coefficient of cos(tau) in the forcing term is zero.
    # This will give us an equation to solve for lmbda2.
    secular_eq_l2 = forcing_R3_expanded.coeff(cos(tau))
    
    # Solve for lmbda2
    lmbda2_sol = sympy.solve(secular_eq_l2, lmbda2)[0]
    lmbda2_sol = sympy.simplify(lmbda2_sol)

    # We have lambda = omega^2. So lmbda = lmbda0 + eps*lmbda1 + eps^2*lmbda2 + ...
    # omega = sqrt(lmbda) = sqrt(lmbda0) * sqrt(1 + eps*lmbda1/lmbda0 + eps^2*lmbda2/lmbda0 + ...)
    # omega = omega0 * (1 + 1/2*(eps*lmbda1/lmbda0 + eps^2*lmbda2/lmbda0) + ...)
    # omega = omega0 + eps*(omega0*lmbda1/(2*lmbda0)) + eps^2*(omega0*lmbda2/(2*lmbda0)) + ...
    # omega = omega0 + eps*omega1 + eps^2*omega2 + ...
    # So, omega2 = omega0 * lmbda2 / (2*lmbda0) = sqrt(lmbda0) * lmbda2 / (2*lmbda0) = lmbda2 / (2*sqrt(lmbda0))
    
    omega0_sol = sqrt(lmbda0_sol)
    omega2_sol = lmbda2_sol / (2 * omega0_sol)
    
    final_expression = sympy.simplify(omega2_sol)
    
    # The question asks for the "3rd term of the nonlinear correction".
    # The correction is eps*omega1 + eps^2*omega2 + eps^3*omega3 + ...
    # Since omega1=0, the first non-zero term is at order eps^2.
    # We interpret the question as asking for omega_2, the coefficient of the first non-zero correction term.
    
    print("The coefficient of the first non-zero nonlinear frequency correction term (omega_2) is:")
    sympy.pprint(final_expression)

solve_rayleigh_plesset_correction()