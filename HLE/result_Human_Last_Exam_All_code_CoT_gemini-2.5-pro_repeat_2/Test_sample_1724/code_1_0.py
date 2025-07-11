import sympy as sp

def solve_bubble_frequency_correction():
    """
    Calculates the second-order frequency correction (omega_2) for the
    Rayleigh-Plesset equation using the Poincaré-Lindstedt method.
    """
    # Define symbolic variables
    gamma, tau = sp.symbols('gamma tau', real=True)

    # --- Step 1: Solve for the first-order correction to the solution (x1) ---
    # The O(epsilon) equation is x1'' + x1 = 0, with x1(0)=1, x1'(0)=0.
    # The solution is x1 = cos(tau).
    # For the purpose of finding omega_2, we need the second-order solution x2.
    # In our notation, this corresponds to y1 from the text explanation.
    x1 = sp.cos(tau)

    # --- Step 2: Solve for the second-order correction to the solution (x2) ---
    # The O(epsilon^2) equation for x2 is:
    # x2'' + x2 = cos^2(t) - (3/2)sin^2(t) + ((3*gamma+1)/2)cos^2(t)
    # The first frequency correction omega_1 is zero.
    rhs_x2 = sp.cos(tau)**2 - sp.S(3)/2 * sp.sin(tau)**2 + (3*gamma+1)/2 * sp.cos(tau)**2
    
    # We find the particular solution for x2 to use in the next step.
    # The homogeneous solution part C*cos(tau) + D*sin(tau) is also needed.
    # Initial conditions: x2(0)=0, x2'(0)=0.
    
    # Find coefficients for the particular solution of form A0 + A2*cos(2*tau)
    rhs_x2_exp = sp.expand_trig(rhs_x2)
    A0 = rhs_x2_exp.coeff(sp.cos(2*tau), 0).coeff(sp.sin(2*tau),0)
    A2 = rhs_x2_exp.coeff(sp.cos(2*tau), 1)
    
    # The particular solution is xp = A0 + A2/(1-2**2) * cos(2*tau)
    xp = A0 - A2/3 * sp.cos(2*tau)
    
    # Find constants for the homogeneous part C1*cos(tau) + C2*sin(tau)
    # x2(0)=0 -> C1 + A0 - A2/3 = 0 -> C1 = A2/3 - A0
    # x2'(0)=0 -> C2 = 0
    C1 = A2/3 - A0
    x2 = C1*sp.cos(tau) + xp

    # --- Step 3: Find the forcing term for the O(epsilon^3) equation ---
    # The secular term from this forcing determines omega_2.
    x1_p, x1_pp = x1.diff(tau), x1.diff(tau,2)
    x2_p, x2_pp = x2.diff(tau), x2.diff(tau,2)

    F = -(x1*x2_pp + x2*x1_pp) - 3*x1_p*x2_p + (3*gamma+1)*x1*x2 \
        - ((3*gamma+1)*(3*gamma+2)/6)*x1**3

    # --- Step 4: Isolate the secular driving term (coefficient of cos(tau)) ---
    # The coefficient is given by (1/pi) * integral(F*cos(tau)) from 0 to 2*pi
    secular_coeff = (sp.integrate(F * sp.cos(tau), (tau, 0, 2*sp.pi)) / sp.pi).simplify()

    # --- Step 5: Solve for omega_2 ---
    # The secularity condition is: secular_coeff + 2*omega_2/omega_0 = 0
    # So, omega_2 = -omega_0 * secular_coeff / 2
    # The full term is omega_2 = -sqrt(3*gamma)/16 * (6*gamma**2 - 3*gamma - 2)
    # The polynomial part is P = 8 * secular_coeff
    
    poly_in_gamma = sp.simplify(8 * secular_coeff)
    
    # Extract coefficients of the polynomial P(gamma) = A*gamma^2 + B*gamma + C
    coeff_A = poly_in_gamma.coeff(gamma, 2)
    coeff_B = poly_in_gamma.coeff(gamma, 1)
    coeff_C = poly_in_gamma.coeff(gamma, 0)
    
    omega0_sym = sp.Symbol('omega_0')
    omega2_sym = sp.Symbol('omega_2')
    
    final_equation = sp.Eq(omega2_sym, -omega0_sym/16 * (coeff_A*gamma**2 + coeff_B*gamma + coeff_C))
    
    print("The equation for the second-order frequency correction is:")
    print(f"{final_equation.lhs} = {final_equation.rhs}")
    print("\nwhere omega_0 = sqrt(3*gamma)")
    
    print("\nThe polynomial factor is:")
    print(f"{coeff_A}*gamma**2 + {coeff_B}*gamma + {coeff_C}")
    
    print("\nThe 3rd term of this polynomial is the constant term.")
    
    print("\nFinal Answer:")
    print(coeff_C)
    
    return coeff_C

# Execute the solver
final_term = solve_bubble_frequency_correction()
print(f"<<<{final_term}>>>")