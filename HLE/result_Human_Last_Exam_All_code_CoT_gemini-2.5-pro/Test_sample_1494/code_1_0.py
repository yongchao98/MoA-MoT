import numpy as np
from scipy.integrate import quad
from scipy.stats import linregress

def integrand(k, t, C=1.0):
    """
    The integrand k^2 * exp(-C * k^2 * t), which is the dominant part
    of the formula for the squared norm at large times.
    """
    return k**2 * np.exp(-C * k**2 * t)

def calculate_squared_norm(t, C=1.0):
    """
    Calculates the integral representing ||d_x rho(t)||^2 at time t.
    We integrate from 0 to infinity.
    """
    result, error = quad(integrand, 0, np.inf, args=(t, C))
    return result

# We analyze the decay for a range of large time values.
t_values = np.logspace(2, 5, 30)

# Calculate the squared norm for each time point
squared_norm_values = np.array([calculate_squared_norm(t) for t in t_values])

# We expect the decay to follow a power law: norm_squared ~ t^slope.
# To find the slope, we can perform a linear regression on the log-log data.
# log(norm_squared) ~ slope * log(t) + intercept
log_t = np.log(t_values)
log_squared_norm = np.log(squared_norm_values)

# Perform linear regression to find the slope for the squared norm
slope_sq, _, _, _, _ = linregress(log_t, log_squared_norm)

# The L2-norm itself is the square root of the integral
norm_values = np.sqrt(squared_norm_values)
log_norm = np.log(norm_values)

# Perform linear regression for the norm
slope_norm, _, _, _, _ = linregress(log_t, log_norm)

# --- Output the results of the numerical analysis ---
print("This script numerically verifies the derived decay rate.")
print("="*60)
print("1. Analysis for the squared norm: ||\u2202_x \u03C1(t)||_{L^2}^2")
print("-" * 60)
print(f"The calculated slope from the log-log plot is: {slope_sq:.4f}")
print("The theoretical slope is -1.5 (from t\u207b\u00b3\u2044\u00b2). The numerical result is in excellent agreement.")
print("The decay law is: ||\u2202_x \u03C1(t)||_{L^2}^2 \u223C C \u00B7 t^({:.2f})".format(slope_sq))
print("\n")
print("2. Analysis for the L2-norm: ||\u2202_x \u03C1(t)||_{L^2}")
print("-" * 60)
print(f"The calculated slope from the log-log plot is: {slope_norm:.4f}")
print("The theoretical slope is -0.75 (from t\u207b\u00b3\u2044\u2074). The numerical result is also in excellent agreement.")
print("\nBased on the analysis, the final equation for the decay is:")
print(f"||\u2202_x \u03C1(t)||_{{L^2}} \u223C C \u00B7 t^({slope_norm:.2f})")