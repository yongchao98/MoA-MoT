import math

# Plan:
# 1. Define the parameters for the Gaussian-based IT3 MF.
# 2. Define the primary UMF and LMF (Gaussian functions) that form the reference Type-2 MF.
# 3. Choose a specific value for the primary input variable 'x' to define the vertical slice.
# 4. Calculate the values of the primary UMF and LMF at this specific 'x'.
# 5. Construct and print the final mathematical formulations for the upper and lower bounds
#    of the vertical cross-section, substituting the calculated values.
#    The formulation will be a function of the secondary variable 'u'.

# --- Step 1: Define Parameters ---
# Parameters for the primary Gaussian membership functions
c = 2.0  # Center of the Gaussians
overline_sigma = 1.0  # Standard deviation for the Upper Membership Function (UMF)
underline_sigma = 0.5 # Standard deviation for the Lower Membership Function (LMF)

# Parameters for the 'blurring' in the third dimension
# These control the uncertainty of the vertical cross-section
overline_alpha = 0.3  # Proportionality constant for the UMF's vertical uncertainty
underline_alpha = 0.4 # Proportionality constant for the LMF's vertical uncertainty

# --- Step 2: Define primary UMF and LMF functions ---
def primary_umf(x, c, sigma):
    """Calculates the primary Upper Membership Function value."""
    return math.exp(-0.5 * ((x - c) / sigma)**2)

def primary_lmf(x, c, sigma):
    """Calculates the primary Lower Membership Function value."""
    return math.exp(-0.5 * ((x - c) / sigma)**2)

# --- Step 3: Choose a specific 'x' for the vertical cross-section ---
x_fixed = 1.5

# --- Step 4: Calculate primary UMF and LMF values at x_fixed ---
overline_mu_A_x = primary_umf(x_fixed, c, overline_sigma)
underline_mu_A_x = primary_lmf(x_fixed, c, underline_sigma)

# --- Step 5: Construct and print the mathematical formulation ---
# The vertical cross-section is a Type-2 fuzzy set whose membership function for a given 'u'
# is an interval [underline_f(u), overline_f(u)]. These are the formulations for those bounds.

print("Mathematical Formulation for the Vertical Cross-Section of an IT3 MF")
print("=====================================================================")
print(f"For a fixed primary input x = {x_fixed}, the vertical cross-section is a Type-2 MF.")
print("Its domain for the secondary variable 'u' is [underline_μ_A(x), overline_μ_A(x)].")
print(f"In this case, u ∈ [{underline_mu_A_x:.4f}, {overline_mu_A_x:.4f}].\n")
print("The membership grade at this cross-section is itself an interval, defined by")
print("upper and lower Gaussian-based functions of 'u'.\n")

# --- Formulation for the Upper Bound of the vertical cross-section ---
print("--- Upper Bound Formulation (as a function of u) ---")
# The general form is: overline_f(u) = exp(-0.5 * ((u - overline_μ_A(x)) / (overline_α * overline_μ_A(x)))**2)
mean_upper = overline_mu_A_x
std_dev_upper = overline_alpha * overline_mu_A_x

print("The equation for the upper bound of the vertical cross-section is:")
print(f"overline_f(u) = exp(-0.5 * ((u - {mean_upper:.4f}) / {std_dev_upper:.4f})**2)")


# --- Formulation for the Lower Bound of the vertical cross-section ---
print("\n--- Lower Bound Formulation (as a function of u) ---")
# The general form is: underline_f(u) = exp(-0.5 * ((u - underline_μ_A(x)) / (underline_α * underline_μ_A(x)))**2)
mean_lower = underline_mu_A_x
std_dev_lower = underline_alpha * underline_mu_A_x

print("The equation for the lower bound of the vertical cross-section is:")
print(f"underline_f(u) = exp(-0.5 * ((u - {mean_lower:.4f}) / {std_dev_lower:.4f})**2)")