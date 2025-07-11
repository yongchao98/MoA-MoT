import sympy as sp

def solve_blowup_condition():
    """
    Analyzes the system of ODEs to find the blow-up condition for y(0).
    """
    # Define variables
    x, y = sp.symbols('x y')
    x0 = sp.Symbol('x0')

    # --- Part 1: Qualitative Analysis ---
    # The system is x'(t) = -3x(t)y(t), y'(t) = -y^2(t) - x(t) + 1.
    # We are given x(0) = x0 > 1.
    #
    # Case 1: y(0) < 0.
    # If y(t) < 0, then x'(t) = -3xy > 0, so x(t) is increasing.
    # Since x(0) > 1, x(t) will remain > 1.
    # y'(t) = -y^2 - (x(t) - 1). Since x(t) > 1, x(t)-1 > 0.
    # So, y'(t) < -y^2. By comparison with the ODE u' = -u^2, whose solutions are u(t) = 1/(t-c),
    # y(t) must go to -infinity in finite time. So, if y(0) < 0, the solution blows up.
    #
    # Case 2: y(0) = 0.
    # At t=0, y'(0) = -0^2 - x0 + 1 = 1 - x0 < 0 (since x0 > 1).
    # This means y(t) becomes negative for t > 0.
    # The trajectory enters the y < 0 region, leading to blow-up as in Case 1.
    #
    # Case 3: y(0) > 0.
    # We need to determine if the trajectory crosses y=0 (blow-up) or x=1.
    # This is determined by the stable manifold of the saddle point (1,0).

    # --- Part 2 & 3: Find Separatrix Equation ---
    # The ODE for trajectories is dy/dx = y'/x' = (-y**2 - x + 1) / (-3*x*y)
    # This can be written as (y**2 + x - 1)dx + (3*x*y)dy = 0.
    # Let's find the potential function Phi(x,y).
    
    M = y**2 + x - 1
    N = 3*x*y

    # Find integrating factor mu(x)
    integrand = (sp.diff(M, y) - sp.diff(N, x)) / N
    mu = sp.exp(sp.integrate(integrand, x)) # mu(x) = x**(-1/3)

    # The exact ODE is (mu*M)dx + (mu*N)dy = 0
    M_new = sp.simplify(M * mu)
    N_new = sp.simplify(N * mu)

    # Integrate to find Phi(x,y)
    Phi = sp.integrate(N_new, y) + sp.integrate(M_new - sp.diff(sp.integrate(N_new, y), x), x)
    Phi = sp.simplify(Phi)
    # Phi = 3*x**(2/3)*y**2/2 + 3*x**(5/3)/5 - 3*x**(2/3)/2

    # The separatrix passes through the saddle point (1, 0).
    C_sep = Phi.subs({x: 1, y: 0})

    # Print the equation of the separatrix
    # We manually format the fractions and exponents for clarity.
    print("The equation for the trajectories is given by Phi(x, y) = C, where:")
    print("Phi(x, y) = (3/2)*x^(2/3)*y^2 + (3/5)*x^(5/3) - (3/2)*x^(2/3)")
    print(f"The separatrix is the trajectory passing through the saddle point (1, 0), so its constant is C = {C_sep}.")
    
    # We use rational numbers for exact representation in the final printed equation.
    n1, d1 = sp.S(3)/2, sp.S(3)/2
    n2, d2 = sp.S(3)/5, sp.S(5)/3
    n3, d3 = sp.S(3)/2, sp.S(2)/3
    C_val = C_sep
    
    print(f"The equation of the separatrix is: ({n1.p}/{n1.q})*x^({d3.p}/{d3.q})*y^2 + ({n2.p}/{n2.q})*x^({d2.p}/{d2.q}) - ({n3.p}/{n3.q})*x^({d3.p}/{d3.q}) = {C_val.p}/{C_val.q}")
    print("-" * 20)

    # --- Part 4 & 5: Analyze Blow-up for y(0) > 0 ---
    # For y(0) > 0, blow-up occurs if Phi(x0, y0) < C_sep.
    # This leads to the inequality y0**2 < f(x0). Let's find f(x).
    
    # y0^2 < (C_sep - (3/5*x0**(5/3) - 3/2*x0**(2/3))) / (3/2*x0**(2/3))
    f_x = sp.simplify((C_sep - (sp.S(3)/5 * x**(sp.S(5)/3) - sp.S(3)/2 * x**(sp.S(2)/3))) / (sp.S(3)/2 * x**(sp.S(2)/3)))
    
    print("For y(0) > 0, blow-up would require the initial condition (x(0), y(0)) to satisfy:")
    # Print the inequality
    print(f"y(0)^2 < 1 - (2/5)*x(0) - (3/5)*x(0)^(-2/3)")
    
    # Analyze the RHS of the inequality for x > 1.
    f_1 = f_x.subs(x, 1)
    df_dx = sp.diff(f_x, x)
    
    print("\nLet f(x) be the right hand side of the inequality.")
    print(f"At x = 1, f(1) = {f_1}.")
    print(f"The derivative is f'(x) = {df_dx}.")
    print("For x > 1, f'(x) is negative, so f(x) < 0.")
    print("Since y(0)^2 must be non-negative, this inequality can never be satisfied for x(0) > 1 and y(0) > 0.")
    print("-" * 20)
    
    print("Conclusion: Blow-up occurs if and only if y(0) is less than or equal to 0.")


if __name__ == '__main__':
    solve_blowup_condition()
    print("\n<<<y(0) <= 0>>>")
