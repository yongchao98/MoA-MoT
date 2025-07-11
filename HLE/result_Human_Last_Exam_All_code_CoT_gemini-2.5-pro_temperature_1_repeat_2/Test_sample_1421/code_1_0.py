import sys

# This script prints the formula for the fermionic partition function.
# It does not perform any calculations.

# --- Define components of the formula ---
partition_function = "Z"
equals = "="
path_integral = "∫ [Dψ̄ Dψ]"
exponent_open = "exp("
minus_sign = "-"
time_integral_symbol = "∫"
lower_time_limit = "0"
upper_time_limit = "β"
time_differential = "dτ"
space_integral = "∫ d^d x"
lagrangian = "ψ̄(x, τ) [∂_τ + H₀] ψ(x, τ)"
exponent_close = ")"

# --- Print the main formula ---
print("The formula for the fermionic partition function (Z) in the imaginary time path integral formalism is:")
print("")
# The instruction is to print each number in the final equation.
# The explicit number in the standard form of the equation is the lower limit of the time integral, which is 0.
final_equation = (
    f"{partition_function} {equals} {path_integral} "
    f"{exponent_open}{minus_sign} {time_integral_symbol}_{{{lower_time_limit}}}^{{{upper_time_limit}}} {time_differential} "
    f"{space_integral} {lagrangian}{exponent_close}"
)
print(final_equation)
print("")

# --- Print the crucial boundary condition ---
print("This integral is taken over all field configurations satisfying the anti-periodic boundary condition in imaginary time (τ):")
boundary_condition = f"ψ(τ = {upper_time_limit}, x) = {minus_sign}ψ(τ = {lower_time_limit}, x)"
print(boundary_condition)
print("")

# --- Print definitions of the symbols ---
print("Where:")
print(f"  - Z: The partition function.")
print(f"  - ∫ [Dψ̄ Dψ]: The functional path integral over all anticommuting Grassmann fields ψ and ψ̄.")
print(f"  - β: The inverse temperature, defined as 1/(k_B T).")
print(f"  - τ: Imaginary time, which is integrated from the number {lower_time_limit} to β.")
print(f"  - H₀: The single-particle Hamiltonian operator (e.g., for a free particle, H₀ = -ħ²∇²/(2m) - μ, where μ is the chemical potential).")
