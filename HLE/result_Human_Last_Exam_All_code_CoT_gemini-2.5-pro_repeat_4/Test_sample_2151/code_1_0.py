import sympy
from sympy import Symbol, exp, tanh, sech, diff, S

def solve_pde_quantity():
    """
    Solves for the quantity -u(0,1)/2 by using a Taylor series expansion
    of u(0,t) around t=0.
    """
    # Define the symbolic variable
    x = Symbol('x')
    t = Symbol('t')

    # Step 1 & 2: Define initial conditions from the problem statement
    u_x_0 = -2 + (1 - tanh(x)) / (exp(x) + 1)
    du_dt_x_0 = S(1)/4 * (tanh(x) - 1) * sech(x/2)**2 * (tanh(x) - sech(x) - 2)

    # Step 3: Calculate u(0,0)
    u_0_0 = u_x_0.subs(x, 0)
    print(f"Calculated u(0,0) = {u_0_0}")

    # Step 4: Calculate du/dt(0,0)
    du_dt_0_0 = du_dt_x_0.subs(x, 0)
    print(f"Calculated du/dt(0,0) = {du_dt_0_0}")

    # Step 5: Calculate spatial derivatives at (0,0)
    # First spatial derivative
    du_dx_x_0 = diff(u_x_0, x)
    du_dx_0_0 = du_dx_x_0.subs(x, 0)
    print(f"Calculated du/dx(0,0) = {du_dx_0_0}")

    # Second spatial derivative
    d2u_dx2_x_0 = diff(du_dx_x_0, x)
    d2u_dx2_0_0 = d2u_dx2_x_0.subs(x, 0)
    print(f"Calculated d^2u/dx^2(0,0) = {d2u_dx2_0_0}")

    # Step 6: Use the PDE to find d^2u/dt^2(0,0)
    # PDE: du/dt + 1/8*d2u/dt2 + u*du/dx - 1/8*d2u/dx2 - (u-1)u(u+2) = 0
    # Rearranging for d2u/dt2:
    # d2u/dt2 = -8*(du/dt + u*du/dx - (u-1)u(u+2)) + d2u/dx2
    
    u_val = u_0_0
    du_dt_val = du_dt_0_0
    du_dx_val = du_dx_0_0
    d2u_dx2_val = d2u_dx2_0_0
    
    term_nonlinear = (u_val - 1) * u_val * (u_val + 2)
    
    d2u_dt2_0_0 = -8 * (du_dt_val + u_val * du_dx_val - term_nonlinear) + d2u_dx2_val
    print(f"Calculated d^2u/dt^2(0,0) = {d2u_dt2_0_0}")

    # Step 7: Construct the Taylor series for u(0,t) up to the t^2 term
    u_0_t = u_0_0 + t * du_dt_0_0 + (t**2 / 2) * d2u_dt2_0_0
    print(f"\nTaylor expansion for u(0,t) is: {u_0_t}")

    # Step 8: Evaluate u(0,1)
    u_0_1 = u_0_t.subs(t, 1)
    print(f"Value of u(0,1) = {u_0_0} + ({du_dt_0_0})*1 + ({d2u_dt2_0_0})*1^2/2 = {u_0_1}")

    # Step 9: Calculate the final quantity
    final_quantity = -u_0_1 / 2
    print(f"\nThe final quantity -u(0,1)/2 = -({u_0_1})/2 = {final_quantity}")
    
    # Final answer in the required format
    # The final equation is -(-1/2)/2 = 1/4
    num1 = u_0_1
    num2 = 2
    result = final_quantity
    print(f"The final equation is -({num1})/{num2} = {result}")


solve_pde_quantity()