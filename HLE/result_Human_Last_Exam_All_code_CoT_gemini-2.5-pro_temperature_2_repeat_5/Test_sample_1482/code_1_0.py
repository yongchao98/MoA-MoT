import math

# Define the parameters for the IT3 MF vertical cross-section calculation
# --- Primary Input Variables ---
x = 6.0  # Fixed primary input value
u = 0.8  # Secondary variable, representing the membership level being evaluated

# --- Reference Upper Membership Function (UMF) Parameters ---
# This is the Gaussian function μ_upper(x) that defines the upper bound of the reference Type-2 set
c_x = 5.0     # The center (mean) of the Gaussian UMF
sigma_x = 2.0 # The standard deviation of the Gaussian UMF
height = 1.0  # The peak value or scale of the Gaussian UMF (typically 1.0)

# --- Vertical Slice Parameters ---
# This controls the "blur" or uncertainty in the vertical direction
k_U = 0.1 # Proportionality constant for the standard deviation of the vertical slice

# --- Explanation of the model ---
print("This script calculates the membership grade for a vertical cross-section of a Gaussian Interval Type-3 Membership Function (IT3 MF).")
print("\nThe mathematical formulation for the upper bound of the vertical cross-section is:")
print("  μ_U(u|x) = exp(-0.5 * ((u - μ_upper(x)) / σ_U(x))^2)\n")
print("This formula tells us the possibility (our result) that the true membership grade is 'u', given the primary input 'x'.")
print("It is dependent on two key components:")
print("  1. μ_upper(x): The reference Upper Membership Function (UMF) from a Type-2 system.")
print("  2. σ_U(x): The standard deviation of the vertical slice, which models the Type-3 uncertainty.")
print("-" * 50)

# --- Step-by-step calculation ---

# Step 1: Calculate the reference UMF value, μ_upper(x), for the given x
# This is the center of the Gaussian vertical slice.
mu_upper_x_val = height * math.exp(-0.5 * ((x - c_x) / sigma_x)**2)

# Step 2: Calculate the standard deviation of the vertical slice, σ_U(x)
# This models the "blur". We model it as proportional to the reference UMF's height at x.
sigma_U_x_val = k_U * mu_upper_x_val

# Step 3: Calculate the final membership grade of the vertical cross-section
exponent_term = -0.5 * ((u - mu_upper_x_val) / sigma_U_x_val)**2
final_result = math.exp(exponent_term)

# --- Display the full equation with numerical substitutions ---

print("CALCULATION STEPS:")
print(f"\nGiven inputs: x = {x}, u = {u}")
print(f"Reference UMF parameters: c_x = {c_x}, σ_x = {sigma_x}, height = {height}")
print(f"Vertical slice uncertainty factor: k_U = {k_U}")

print("\n1. Calculate the reference UMF value, μ_upper(x):")
print(f"   μ_upper(x={x}) = {height} * exp(-0.5 * (({x} - {c_x}) / {sigma_x})^2) = {mu_upper_x_val:.6f}")

print("\n2. Calculate the vertical slice standard deviation, σ_U(x):")
print(f"   σ_U(x={x}) = {k_U} * μ_upper(x) = {k_U} * {mu_upper_x_val:.6f} = {sigma_U_x_val:.6f}")

print("\n3. Calculate the final result using the full formulation, μ_U(u|x):")
print(f"\n   μ_U(u={u} | x={x}) = exp(-0.5 * (({u} - {mu_upper_x_val:.6f}) / {sigma_U_x_val:.6f})^2)")
print(f"                  = exp({exponent_term:.6f})")
print(f"                  = {final_result:.6f}")

# The final answer in the requested format
final_equation = f"μ_U(u={u}|x={x}) = exp(-0.5 * (({u} - {mu_upper_x_val:.4f}) / {sigma_U_x_val:.4f})^2) = {final_result:.4f}"
print(f"\n<<<{final_equation}>>>")
