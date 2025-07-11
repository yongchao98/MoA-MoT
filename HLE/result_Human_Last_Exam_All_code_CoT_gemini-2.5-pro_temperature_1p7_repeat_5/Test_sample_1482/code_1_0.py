import math

# --- Plan Explanation ---
# This script determines the mathematical formulation for a vertical cross-section
# of a Gaussian-based Interval Type-3 Membership Function (IT3 MF).
# The steps are detailed in the print statements below.

print("Deriving the formulation for a vertical cross-section of a Gaussian IT3 MF.\n")

# Step 1: Define parameters for the Gaussian IT3 MF model
c = 5.0
sigma_L = 1.0  # Standard deviation for the Lower Membership Function (LMF)
sigma_U = 1.5  # Standard deviation for the Upper Membership Function (UMF)
x_val = 6.0    # The fixed value of the primary variable 'x' for the slice
N = 6.0        # A scaling constant for the standard deviation of the vertical slice

print("--- Step 1: Define Primary UMF and LMF ---")
print("The UMF and LMF for the primary variable 'x' are defined as Gaussian functions.")
print(f"A common model uses a fixed center 'c' and uncertain standard deviation in [sigma_L, sigma_U].")
print(f"UMF: overline_mu_A(x) = exp(-0.5 * ((x - c) / sigma_U)^2)")
print(f"LMF: underline_mu_A(x) = exp(-0.5 * ((x - c) / sigma_L)^2)")
print(f"Model parameters: c = {c}, sigma_L = {sigma_L}, sigma_U = {sigma_U}")
print("-" * 60)

print(f"--- Step 2: Calculate UMF and LMF values at x = {x_val} ---")
# Step 2: Calculate mu_upper and mu_lower at x_val
mu_upper = math.exp(-0.5 * ((x_val - c) / sigma_U)**2)
mu_lower = math.exp(-0.5 * ((x_val - c) / sigma_L)**2)

print(f"Evaluating at x = {x_val}:")
print(f"  overline_mu_A({x_val}) = exp(-0.5 * (({x_val} - {c}) / {sigma_U})^2) = {mu_upper:.4f}")
print(f"  underline_mu_A({x_val}) = exp(-0.5 * (({x_val} - {c}) / {sigma_L})^2) = {mu_lower:.4f}")
print(f"\nThis defines the uncertain membership domain for the secondary variable 'u' as the interval: [{mu_lower:.4f}, {mu_upper:.4f}]")
print("-" * 60)

print("--- Step 3: Formulate the Vertical Slice's Secondary Membership Function ---")
print("The vertical slice at x' is a Type-2 MF. Its upper bound, overline_mu_{A,x'}(u),")
print("is modeled as a Gaussian function of the secondary variable 'u'.")
print("Its general form is: exp(-0.5 * ((u - center_u) / std_dev_u)^2)")

# Step 3: Calculate parameters for the vertical slice Gaussian
center_u = (mu_upper + mu_lower) / 2.0
std_dev_u = (mu_upper - mu_lower) / N

print("\nCalculating the parameters for this Gaussian:")
print(f"  center_u = (overline_mu_A + underline_mu_A) / 2 = ({mu_upper:.4f} + {mu_lower:.4f}) / 2 = {center_u:.4f}")
print(f"  std_dev_u = (overline_mu_A - underline_mu_A) / N = ({mu_upper:.4f} - {mu_lower:.4f}) / {N} = {std_dev_u:.4f}")
print("-" * 60)

print("--- Step 4: Final Formulation ---")
print("Substituting the calculated values, the mathematical formulation for the upper bound")
print(f"of the vertical cross-section at x = {x_val} is:")

# Step 4: Construct and print the final equation string
final_equation = f"overline_mu_{{A,x={x_val}}}(u) = exp(-0.5 * ((u - {center_u:.4f}) / {std_dev_u:.4f})^2)"
print(f"\n  {final_equation}\n")