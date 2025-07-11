import numpy as np
from scipy.special import lambertw

def calculate_total_charge():
    """
    Calculates the total charge on a perturbed spherical droplet.

    The problem asks to calculate the total charge Q by integrating the surface
    charge density sigma over the surface area of a perturbed sphere.
    
    The formula for the total charge Q is given by the integral:
    Q = integral(integral(sigma(theta, phi) * dA))
    
    Substituting the given expressions for sigma and dA, we have:
    sigma * dA = (sigma0 * theta / sin(theta)) * (1/R) * (W_term) * (R^2 * sin(theta) d(theta) d(phi))
    where W_term = W(exp(q_i*theta*phi)) / (1 + W(exp(q_i*theta*phi)))**3.
    
    The terms R and sin(theta) cancel out, simplifying the integrand to:
    sigma * dA = sigma0 * theta * R * W_term * d(theta) * d(phi)
    
    The radius R is R = R0 * (1 + epsilon * sin(n*theta) * cos(m*phi)). Since the
    parameters epsilon, n, and m are not provided, we assume that the integral
    of the perturbation term evaluates to zero. This is a standard approach
    when unspecified parameters are present in such problems. Thus, we can replace R
    with its average value R0 for the calculation.
    
    The integral becomes:
    Q = sigma0 * R0 * integral from 0 to 2*pi d(phi) [ integral from 0 to pi d(theta) (theta * W_term) ]
    
    Using the identity that the term (theta * W_term) is the partial derivative with 
    respect to phi of a simpler function:
    theta * W_term = -1/q_i * d/d(phi) [1 / (1 + W(exp(q_i*theta*phi)))]
    
    We integrate with respect to phi first, which simplifies to:
    Integral_phi = -1/q_i * [1/(1+W(exp(q_i*theta*2*pi))) - 1/(1+W(exp(0)))]
               = 1/q_i * [1/(1+W(1)) - 1/(1+W(exp(2*pi*q_i*theta)))]
               
    Then we integrate the above expression with respect to theta from 0 to pi.
    The integral of the second term is found using the identity:
    integral( d(x) / (1 + W(exp(x))) ) = ln(W(exp(x)))
    
    This leads to the final analytical formula for Q.
    """
    # Given constants
    sigma0 = 7.43e-7  # units of e/nm
    R0 = 30.0         # units of nm
    qi = 2 * np.pi    # dimensionless

    # Physical constants and helper values
    pi = np.pi
    omega = lambertw(1).real # W(1) is a real number
    
    # Final derived formula for Q
    # Q = (sigma0 * R0 / qi) * [pi / (1 + omega) - (ln(W(exp(2*pi^2*qi))) - ln(omega)) / (2*pi*qi)]
    
    print("Step 1: Define the constants")
    print(f"sigma0 = {sigma0:.2e} e/nm")
    print(f"R0 = {R0:.1f} nm")
    print(f"qi = 2 * pi = {qi:.4f}")
    print(f"omega = W(1) = {omega:.4f}")
    print("-" * 30)

    # --- Calculation breakdown ---
    print("Step 2: Calculate each term in the formula for Q")
    # Formula: Q = A * (B - C)
    # A = sigma0 * R0 / qi
    # B = pi / (1 + omega)
    # C = D / E where D = ln(W(exp(2*pi^2*qi))) - ln(omega) and E = 2*pi*qi

    A = sigma0 * R0 / qi
    print(f"Term A (pre-factor) = (sigma0 * R0 / qi) = ({sigma0:.2e} * {R0:.1f}) / {qi:.4f} = {A:.4e}")

    B = pi / (1 + omega)
    print(f"Term B (first term in brackets) = pi / (1 + omega) = {pi:.4f} / (1 + {omega:.4f}) = {B:.4f}")
    
    exp_arg = 2 * pi**2 * qi
    W_large_arg = np.exp(exp_arg)
    # Use real part of Lambert W for real inputs
    W_large_val = lambertw(W_large_arg).real 
    
    log_W_large = np.log(W_large_val)
    log_omega = np.log(omega)
    
    D = log_W_large - log_omega
    E = 2 * pi * qi
    C = D / E
    
    print(f"Calculating Term C:")
    print(f"  Exponent for W function = 2*pi^2*qi = {exp_arg:.4f}")
    print(f"  W(exp({exp_arg:.4f})) = {W_large_val:.4f}")
    print(f"  Numerator D = ln({W_large_val:.4f}) - ln({omega:.4f}) = {D:.4f}")
    print(f"  Denominator E = 2*pi*qi = {E:.4f}")
    print(f"Term C (second term in brackets) = D / E = {C:.4f}")
    print("-" * 30)

    print("Step 3: Combine terms to find the total charge Q")
    
    total_charge = A * (B - C)
    
    # Final equation printout
    print("Q = (sigma0 * R0 / qi) * [pi / (1 + omega) - (ln(W(exp(2*pi^2*qi))) - ln(omega)) / (2*pi*qi)]")
    print(f"Q = {A:.4e} * [{B:.4f} - {C:.4f}]")
    print(f"Q = {A:.4e} * [{B - C:.4f}]")
    print(f"Q = {total_charge:.4f} e")

    return total_charge

# Execute the calculation
final_charge = calculate_total_charge()
print("\nFinal Answer:")
print(f"The total charge on the droplet is {final_charge:.2f} e.")
# A specific format for the answer is requested.
print(f"<<<{final_charge:.2f}>>>")
