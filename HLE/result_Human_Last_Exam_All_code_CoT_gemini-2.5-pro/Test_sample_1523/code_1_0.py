# This script prints the final expression for the inner product (φ, D_φ).

# The numerical coefficients in the final, simplified equation are all 1.
coeff1 = 1
coeff2 = 1
coeff3 = 1

# Construct the output string for the equation.
# We use unicode characters for a better mathematical representation.
final_equation = f"(φ, Dφ) = ∫_0^β dτ ∫ dᵈx⃗  [ {coeff1}⋅(∂_τ φ)² + {coeff2}⋅(∇φ)² + {coeff3}⋅m²φ² ]"

print("The inner product (φ, D_φ) for a neutral scalar field is given by the following equation:")
print("-" * 70)
print(final_equation)
print("-" * 70)
print("\nWhere:")
print("  - φ = φ(τ, x⃗) is the scalar field.")
print("  - τ is the imaginary time, integrated from 0 to β (inverse temperature).")
print("  - x⃗ represents the d spatial dimensions, integrated over all space.")
print("  - (∂_τ φ)² is the square of the derivative with respect to imaginary time.")
print("  - (∇φ)² is the square of the spatial gradient.")
print("  - m is the mass of the field.")