import math

def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the first-order epsilon expansion.
    """
    # Spatial dimension
    d = 3
    
    # Upper critical dimension for the n-vector model
    d_c = 4
    
    # Number of components of the order parameter (n=1 for Ising model)
    n = 1
    
    # Calculate epsilon (ϵ)
    epsilon = d_c - d
    
    # Calculate the terms for the first-order alpha formula
    # Formula: α = (4 - n) * ϵ / (2 * (n + 8))
    numerator = 4 - n
    denominator_term_1 = 2
    denominator_term_2 = n + 8
    
    # Calculate the final value of alpha
    alpha = (numerator * epsilon) / (denominator_term_1 * denominator_term_2)
    
    # --- Outputting the explanation and calculation ---
    print("Calculating the specific heat exponent α for d=3 using the first-order epsilon expansion.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The spatial dimension is d = {d}.")
    print(f"Therefore, epsilon is ϵ = d_c - d = {d_c} - {d} = {epsilon}.")
    print(f"The universality class is assumed to be the Ising model, so n = {n}.")
    
    print("\nThe first-order expansion formula for α is: α = (4 - n) * ϵ / (2 * (n + 8))")
    
    print("\nSubstituting the values into the equation:")
    print(f"α = (4 - {n}) * {epsilon} / (2 * ({n} + 8))")
    print(f"α = ({numerator}) * {epsilon} / ({denominator_term_1} * {denominator_term_2})")
    print(f"α = {numerator * epsilon} / {denominator_term_1 * denominator_term_2}")
    print(f"α = 1/6")
    print(f"The numerical value is α ≈ {alpha}")

if __name__ == "__main__":
    calculate_alpha_exponent()
