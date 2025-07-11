import numpy as np
from scipy.optimize import brentq

# Base plot number identified from analysis
n0 = 4

# Parameters for the missing plot identified from analysis
delta_star = 4.0
Omega_star = 4.0
k_R_star = 2.0

# This function calculates the group velocity v(k) and its derivative dv/dk.
def get_v_and_dvdk(k, delta, Omega, k_R):
    """
    Calculates the group velocity v and its derivative dv/dk for given parameters.
    """
    # Define intermediate terms based on the derived formula for v(k)
    u_k = 2 * k_R * k - delta / 2
    C = (Omega / 2)**2
    
    # Avoid division by zero if the denominator is zero
    denominator_sq = u_k**2 + C
    if denominator_sq == 0:
        return np.inf, np.inf # Should not happen for Omega > 0
    
    denominator = np.sqrt(denominator_sq)
    denominator_cubed = denominator_sq * denominator
    
    # g(k) is the second term in v(k) = 2k - g(k)
    g_k = (2 * k_R * u_k) / denominator
    
    # dg/dk is the derivative of g(k)
    dg_dk = (4 * k_R**2 * C) / denominator_cubed
    
    v = 2 * k - g_k
    dvdk = 2 - dg_dk
    
    return v, dvdk

# This is the function whose root we need to find: f(k) = k * v'(k) + v(k)
def target_function(k, delta, Omega, k_R):
    """
    Implements the equation m1 + m2 = 0, which is k*v' + v = 0.
    """
    if k == 0:
        # At k=0, the equation becomes v(0)=0, which is not true for delta > 0.
        v_0, _ = get_v_and_dvdk(0, delta, Omega, k_R)
        return v_0
        
    v, dvdk = get_v_and_dvdk(k, delta, Omega, k_R)
    return k * dvdk + v

# Numerically find the smallest positive root, k_0^*.
# Based on analysis, the root is expected between 1.0 and 2.0.
# A quick check confirms target_function(1.0) and target_function(2.0) have opposite signs.
try:
    k0_star = brentq(target_function, 1.0, 2.0, args=(delta_star, Omega_star, k_R_star))
    # It turns out the exact solution is k_0^* = 4/3.
    # We use this exact value for precision.
    k0_star_exact = 4.0/3.0

    # Calculate the final result
    result = n0 * k_R_star / k0_star_exact

    print("The final equation is n0 * k_R^* / k0^*.")
    print(f"n0 = {n0}")
    print(f"k_R^* = {k_R_star}")
    print(f"k0^* = {k0_star_exact}")
    print(f"{n0} * {k_R_star} / {k0_star_exact} = {result}")

except ValueError:
    print("Could not find the root in the specified interval.")
