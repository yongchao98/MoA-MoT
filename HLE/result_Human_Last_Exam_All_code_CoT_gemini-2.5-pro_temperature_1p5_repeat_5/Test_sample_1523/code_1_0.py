import math

def display_scalar_inner_product():
    """
    This script explains and displays the inner product (ϕ, D_ϕ)
    for a neutral scalar field in finite-temperature field theory.
    """
    print("The inner product (ϕ, D_ϕ) is a key component of the action S[ϕ] for a free scalar field.")
    print("The action is S[ϕ] = (1/2) * (ϕ, D_ϕ), and it appears in the partition function Z = ∫ Dϕ exp(-S[ϕ]).\n")
    
    print("-" * 60)
    print("1. Position Space Representation")
    print("-" * 60)
    
    print("In Euclidean spacetime (with imaginary time τ), the inner product is an integral:")
    print("   (ϕ, D_ϕ) = ∫ dτ d³x  ϕ(x) * D * ϕ(x)\n")
    
    print("The operator D is given by:")
    print("   D = [-∂_τ² - ∇² + m²]\n")
    
    print("Let's break down the 'numbers' in this operator expression:")
    print(f"   - Coefficient of the time derivative term (∂_τ²): {-1}")
    print(f"   - Power of the time derivative (∂_τ): {2}")
    print(f"   - Coefficient of the spatial derivative term (∇²): {-1}")
    print(f"   - Power of the spatial derivative (∇): {2}")
    print(f"   - Coefficient of the mass term (m²): {1}")
    print(f"   - Power of the mass (m): {2}\n")

    print("-" * 60)
    print("2. Momentum Space Representation")
    print("-" * 60)

    print("In momentum space, the operator D becomes an algebraic expression dependent on the")
    print("Matsubara frequency (ω_n) and spatial momentum (p).\n")

    print("The Matsubara frequency for a boson is:")
    print("   ω_n = (2 * π * n) / β   (where β=1/T, n is an integer)\n")

    print("The 'numbers' in the Matsubara frequency formula are:")
    print(f"   - Integer factor: {2}")
    print(f"   - Mathematical constant: π ≈ {math.pi}\n")

    print("The operator D in momentum space is:")
    print("   D(ω_n, p) = ω_n² + p² + m²\n")

    print("The 'numbers' in this expression are the powers of the variables:")
    print(f"   - Power of the Matsubara frequency (ω_n): {2}")
    print(f"   - Power of the spatial momentum (p): {2}")
    print(f"   - Power of the mass (m): {2}")

if __name__ == "__main__":
    display_scalar_inner_product()