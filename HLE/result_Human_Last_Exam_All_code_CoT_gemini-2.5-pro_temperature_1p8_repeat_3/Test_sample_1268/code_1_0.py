import math

print("This script determines the upper bound for the norm of an ideal in a quadratic number field in relation to the field's covolume.")
print("We interpret the 'maximum norm' (k_k,∞) as the norm N(a) of an ideal 'a'.")
print("This norm is bounded by Minkowski's constant, M_K. Our goal is to relate M_K to the covolume V.")
print("-" * 60)

print("1. Key Formulas for a Quadratic Field (degree n=2)")
print("   - Let K be a quadratic number field, and V be its covolume.")
print("   - Let r2 be the number of pairs of complex embeddings (r2=0 for real, r2=1 for imaginary fields).")
print("   - Minkowski's Bound: M_K = (4/π)^r2 * (2!/2^2) * sqrt(|Δ_K|)")
print("   - Covolume:          V = (1/2^r2) * sqrt(|Δ_K|)")
print("-" * 60)

print("2. Deriving the Relationship")
print("First, we express the discriminant term using the covolume V:")
print("   From V = (1/2^r2) * sqrt(|Δ_K|), we get sqrt(|Δ_K|) = (2^r2) * V.")
print("\nNext, we substitute this into the formula for M_K:")
print("   M_K = (4/π)^r2 * (1/2) * [ (2^r2) * V ]")
print("   M_K = (1/2) * [ (4*2)/π ]^r2 * V")
print("   M_K = (1/2) * (8/π)^r2 * V")
print("-" * 60)

print("3. Analyzing the Two Cases for Quadratic Fields")

# Case 1: Real Quadratic Field (e.g., Q(sqrt(2)))
r2_real = 0
const_real = (1/2) * ((8/math.pi)**r2_real)
print("\nCase A: Real Quadratic Fields (r2 = 0)")
print(f"   The constant C = (1/2) * (8/π)^0 = 1/2 = {const_real:.4f}")
print("   The relationship is: M_K = (1/2) * V")

# Case 2: Imaginary Quadratic Field (e.g., Q(sqrt(-1)))
r2_imag = 1
const_imaginary = (1/2) * ((8/math.pi)**r2_imag)
print("\nCase B: Imaginary Quadratic Fields (r2 = 1)")
print(f"   The constant C = (1/2) * (8/π)^1 = 4/π ≈ {const_imaginary:.4f}")
print("   The relationship is: M_K = (4/π) * V")
print("-" * 60)

print("4. Conclusion: The General Upper Bound")
print("To establish a single upper bound that holds for all quadratic fields, we must choose the larger of the two constants.")
print(f"Comparing the constants C_real = {const_real:.4f} and C_imaginary = {const_imaginary:.4f}, the larger value is 4/π.")
print("\nThus, the general upper bound for the norm k_k,∞ in relation to the covolume V is given by the inequality:")

# Final Equation Output
final_constant_symbolic = "4/π"
final_constant_numeric = 4 / math.pi
print(f"\n   k_k,∞ <= ({final_constant_symbolic}) * V")
print(f"   k_k,∞ <= {final_constant_numeric:.6f} * V")