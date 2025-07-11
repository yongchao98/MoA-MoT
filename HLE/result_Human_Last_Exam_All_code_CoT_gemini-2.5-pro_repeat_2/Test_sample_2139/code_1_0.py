import numpy as np
from scipy.special import gamma
from scipy.integrate import solve_ivp

# Step 1: Define the initial condition y(0) from the problem statement.
# y(0) = (128 * 3^(1/6) * Gamma(2/3))^(-1)
y0 = (128 * np.power(3, 1/6) * gamma(2/3))**(-1)

# Step 2: Define the ODE as a system of first-order ODEs.
# Let y_vec = [y, y'] = [z1, z2]
# Then y_vec' = [y', y''] = [z2, f(t, z1, z2)]
def ode_system(t, y_vec):
    """
    Defines the system of first-order ODEs.
    y_vec[0] is y(t), and y_vec[1] is y'(t).
    The function returns [y'(t), y''(t)].
    """
    y, dy = y_vec
    
    # Avoid numerical issues at the exact singularity, though we start at t_start > 0.
    if t == 0:
        # Based on y'(0)=0. y''(0) is indeterminate but not needed.
        return [0, 0] 
        
    sin_t = np.sin(t)
    cos_t = np.cos(t)
    
    # Check for domain issues of tan(t) and sec(t)
    if cos_t == 0:
        return [np.nan, np.nan]
        
    tan_t = sin_t / cos_t
    sec_t = 1 / cos_t
    
    t4p1 = t**4 + 1
    
    # Define the coefficients of the original ODE: P2*y'' + P1*y' + P0*y = RHS
    P2_t = 4 * t4p1 * tan_t * sec_t
    
    P1_t = 2 * (t4p1 + 2 * tan_t * (t4p1 * tan_t + 8 * t**3) + 1) * sec_t
    
    P0_t = 8 * t**2 * (t + 2 * tan_t * (t * tan_t + 3)) * sec_t
    
    # Check for domain issues of sqrt(sin(t))
    if sin_t < 0:
      RHS_t = np.nan
    else:
      RHS_t = t4p1 * y * np.sqrt(sin_t)
    
    # Calculate y'' from the ODE: y'' = (RHS - P1*y' - P0*y) / P2
    ddy = (RHS_t - P1_t * dy - P0_t * y) / P2_t
    
    return [dy, ddy]

# Step 3: Set up initial conditions for the numerical solver near t=0.
# We start at a small t > 0 to avoid the singularity.
t_start = 1e-8

# Use the asymptotic solution near t=0: y(t) ~ y0 * (1 + t^(3/2)/9)
# From this, we derive y'(t) ~ y0 * (1/6) * t^(1/2)
y_start = y0 * (1 + t_start**1.5 / 9.0)
dy_start = y0 / 6.0 * t_start**0.5

initial_conditions = [y_start, dy_start]

# Step 4: Solve the ODE numerically from t_start to t=pi/4.
t_final = np.pi/4
t_span = [t_start, t_final]

solution = solve_ivp(
    ode_system, 
    t_span, 
    initial_conditions, 
    method='RK45', 
    dense_output=True,
    rtol=1e-9, 
    atol=1e-12
)

# Step 5: Extract and print the final result.
y_pi_over_4 = solution.sol(t_final)[0]

print(f"The initial radius at t=0 is y(0) = {y0}")
print(f"The radius of the balloon at t=pi/4 is y(pi/4) = {y_pi_over_4}")