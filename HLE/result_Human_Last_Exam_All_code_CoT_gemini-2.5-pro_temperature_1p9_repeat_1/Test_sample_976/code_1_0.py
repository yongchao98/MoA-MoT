def print_solution():
    """
    Prints the final expressions for the electric potential and electric field
    for a conductive sphere in a uniform electric field.
    """
    
    potential_expression = r"""
Φ(r, θ) = 
⎧ 
⎪ - (3 * σ₂ / (σ₁ + 2*σ₂)) * E₀ * r * cos(θ)               for r < R
⎨ 
⎪ -E₀ * ( r - ((σ₁ - σ₂) * R³) / ((σ₁ + 2*σ₂) * r²) ) * cos(θ)    for r > R
⎩
"""

    field_expression = r"""
E⃗(r, θ) = 
⎧ 
⎪  (3*σ₂*E₀ / (2*σ₂ + σ₁)) * (cos(θ)r̂ - sin(θ)θ̂)                                                                   for r < R
⎪
⎨ E₀ * [ 1 + (2*(σ₁ - σ₂) * R³) / ((2*σ₂ + σ₁) * r³) ] * cos(θ)r̂ 
⎪  - E₀ * [ 1 - ((σ₁ - σ₂) * R³) / ((2*σ₂ + σ₁) * r³) ] * sin(θ)θ̂      for r > R
⎩
"""

    print("The correct expressions for the electric potential and electric field are:")
    print("-" * 70)
    print("Electric Potential Φ(r, θ):")
    print(potential_expression)
    print("-" * 70)
    print("Electric Field E(r, θ):")
    print(field_expression)
    print("-" * 70)
    print("This corresponds to Choice B.")

# Execute the function to print the solution
print_solution()