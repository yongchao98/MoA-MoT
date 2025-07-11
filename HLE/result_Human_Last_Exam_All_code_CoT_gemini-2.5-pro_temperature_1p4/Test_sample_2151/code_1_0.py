import sympy
import math

# Announce the plan
print("This script solves the problem by using a Taylor series expansion for u(0,t) around t=0.")
print("The necessary derivatives are calculated using the initial conditions and the PDE itself.")
print("-" * 30)

# Step 1: Define symbols and initial conditions
x = sympy.Symbol('x')
t = sympy.Symbol('t')

# Given initial conditions as sympy expressions
u_x0_expr = -2 + (1 - sympy.tanh(x)) / (sympy.exp(x) + 1)
ut_x0_expr = sympy.Rational(1, 4) * (sympy.tanh(x) - 1) * (sympy.sech(x/2)**2) * (sympy.tanh(x) - sympy.sech(x) - 2)

# Evaluate initial conditions at x=0
u_00 = u_x0_expr.subs(x, 0)
ut_00 = ut_x0_expr.subs(x, 0)
print(f"Value of u at (0,0): {u_00}")
print(f"Value of u_t at (0,0): {ut_00}")
print("-" * 30)

# Step 2: Calculate spatial derivatives at (0,0)
ux_x0_expr = sympy.diff(u_x0_expr, x)
uxx_x0_expr = sympy.diff(ux_x0_expr, x)

ux_00 = ux_x0_expr.subs(x, 0)
uxx_00 = uxx_x0_expr.subs(x, 0)
print(f"Value of u_x at (0,0): {ux_00}")
print(f"Value of u_xx at (0,0): {uxx_00}")
print("-" * 30)

# Step 3: Calculate u_tt at (0,0) using the PDE
# PDE: u_t + 1/8*u_tt + u*u_x - 1/8*u_xx - (u-1)*u*(u+2) = 0
# u_tt = -8*u_t - 8*u*u_x + u_xx + 8*(u-1)*u*(u+2)
potential_term = (u_00 - 1) * u_00 * (u_00 + 2)
utt_00 = -8 * ut_00 - 8 * u_00 * ux_00 + uxx_00 + 8 * potential_term

print(f"Calculating u_tt(0,0) from the PDE:")
print(f"u_tt(0,0) = -8*({ut_00}) - 8*({u_00})*({ux_00}) + ({uxx_00}) + 8*(({u_00})-1)*({u_00})*(({u_00})+2)")
print(f"u_tt(0,0) = {-8*ut_00} - {-8*u_00*ux_00} + {uxx_00} + {8*potential_term}")
print(f"Value of u_tt at (0,0): {utt_00}")
print("-" * 30)

# Step 4: Use Taylor series to approximate u(0,1)
t_val = 1
u_01_approx = u_00 + ut_00 * t_val + sympy.Rational(1, 2) * utt_00 * t_val**2
print("Approximating u(0,1) using Taylor expansion: u(0,1) ≈ u(0,0) + u_t(0,0)*t + 0.5*u_tt(0,0)*t^2")
print(f"u(0,1) ≈ {u_00} + {ut_00}*({t_val}) + 0.5*({utt_00})*({t_val})^2")
print(f"Value of u(0,1) ≈ {u_01_approx}")
print("-" * 30)

# Final result calculation
result = -u_01_approx / 2
print("Final Calculation:")
print(f"The required quantity is -u(0,1)/2.")
print(f"Result = -({u_01_approx}) / 2 = {result}")

print("<<<{}>>>".format(float(result)))