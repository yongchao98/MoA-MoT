import sympy as sp
from sympy import tanh, sech, exp

def solve_pde_quantity():
    """
    Calculates the required quantity by determining the coefficients for the
    Taylor series expansion of u(0,t) around t=0.
    """
    # Define the symbolic variable
    x = sp.Symbol('x')

    # --- Step 1 & 2: Define initial conditions ---
    u_x0_expr = -2 + (1 - tanh(x)) / (exp(x) + 1)
    ut_x0_expr = sp.Rational(1, 4) * (tanh(x) - 1) * (sech(x/2)**2) * (tanh(x) - sech(x) - 2)

    # Evaluate ICs at x=0
    u_00 = u_x0_expr.subs(x, 0)
    ut_00 = ut_x0_expr.subs(x, 0)
    
    print(f"Step 1: u(0,0) = {sp.latex(u_x0_expr.subs(x,0))} = {u_00}")
    print(f"Step 2: du/dt(0,0) = {sp.latex(ut_x0_expr.subs(x,0))} = {ut_00}")

    # --- Step 3: Calculate spatial derivatives at (0,0) ---
    ux_x0_expr = sp.diff(u_x0_expr, x)
    uxx_x0_expr = sp.diff(ux_x0_expr, x)

    ux_00 = ux_x0_expr.subs(x, 0)
    uxx_00 = uxx_x0_expr.subs(x, 0)
    
    print(f"Step 3a: du/dx(0,0) = {sp.latex(ux_00)} = {ux_00}")
    print(f"Step 3b: d^2u/dx^2(0,0) = {sp.latex(uxx_00)} = {uxx_00}")

    # --- Step 4: Use the PDE to find d^2u/dt^2(0,0) ---
    # PDE: u_t + 1/8*u_tt + u*u_x - 1/8*u_xx - (u-1)*u*(u+2) = 0
    # Rearranging for u_tt:
    # u_tt = -8 * (u_t + u*u_x - 1/8*u_xx - (u-1)*u*(u+2))
    
    u, ut, ux, uxx = sp.symbols('u ut ux uxx')
    
    # Substitute the calculated numerical values
    term_ut = ut_00
    term_u_ux = u_00 * ux_00
    term_uxx = sp.Rational(1, 8) * uxx_00
    term_nonlinear = (u_00 - 1) * u_00 * (u_00 + 2)

    utt_00 = -8 * (term_ut + term_u_ux - term_uxx - term_nonlinear)
    
    print(f"Step 4: Using the PDE, d^2u/dt^2(0,0) is calculated.")
    print(f"   u_t = {term_ut}")
    print(f"   u*u_x = ({u_00}) * ({ux_00}) = {term_u_ux}")
    print(f"   -1/8*u_xx = -1/8 * ({uxx_00}) = {-term_uxx}")
    print(f"   -(u-1)u(u+2) = -({u_00-1})*({u_00})*({u_00+2}) = {-term_nonlinear}")
    print(f"   d^2u/dt^2(0,0) = {sp.latex(utt_00)} = {utt_00}")
    
    # --- Step 5: Construct and evaluate the Taylor Polynomial for u(0,t) ---
    t = sp.Symbol('t')
    u_0t_poly = u_00 + ut_00 * t + sp.Rational(1, 2) * utt_00 * t**2

    # Evaluate at t=1
    u_01 = u_0t_poly.subs(t, 1)
    
    print("\nStep 5: The Taylor series for u(0,t) is:")
    print(f"   u(0,t) = {sp.latex(u_00)} + ({sp.latex(ut_00)})t + \\frac{{{sp.latex(utt_00)}}}{{2}}t^2 + ...")
    print(f"   u(0,t) = {u_00} + {ut_00}*t + {utt_00/2}*t^2 + ...")
    print(f"   Evaluating at t=1: u(0,1) = {sp.latex(u_00)} + {sp.latex(ut_00)} + {sp.latex(utt_00/2)} = {u_01}")

    # --- Step 6: Compute the final quantity ---
    final_quantity = -u_01 / 2
    
    print("\nStep 6: The final required quantity is -u(0,1)/2:")
    print(f"   -u(0,1)/2 = -({u_01})/2 = {final_quantity}")
    
    print("\nFinal Answer Equation:")
    print(f"u(0,1) = {u_00} + {ut_00}*(1) + {utt_00/2}*(1)^2 = {u_01}")
    print(f"Result = -({u_01}) / 2 = {final_quantity}")

solve_pde_quantity()