import math

def print_force_equation():
    """
    This function prints the derived formula for the instantaneous force f_x(t).
    The formula corresponds to answer choice B.
    """
    
    # Define variables as strings for printing the formula
    two = "2"
    pi = "pi"
    R = "R"
    N = "N"
    mu0 = "mu_0"
    alpha_T = "alpha_T"
    T = "T"
    T0 = "T_0"
    N0 = "N_0"
    I0 = "I_0"
    i0 = "i_0"
    sin_omega_t = "sin(omega * t)"
    g = "g"
    Bs = "B_s"
    
    # Construct the numerator part of the equation
    numerator = f"{two} * {pi} * {R} * {N} * {mu0} * (1 - {alpha_T} * ({T} - {T0})) * {N0} * {I0} * {i0} * {sin_omega_t}"
    
    # Construct the denominator part of the equation
    denominator = f"{g}^2 * (1 + ({mu0} * {N0} * {I0}) / ({g} * {Bs}))"
    
    # Print the final equation with a negative sign
    print(f"f_x(t) = -({numerator}) / ({denominator})")

print_force_equation()