def calculate_scaling_exponent_alpha():
    """
    Calculates the specific heat scaling exponent α using the first-order
    epsilon expansion for a system in d=3 dimensions.
    """
    # Parameters for the calculation
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the n-vector model
    n = 1  # Number of spin components (assuming the Ising model)

    # Epsilon (ϵ) is the expansion parameter
    epsilon = d_c - d

    # Numerator and denominator of the formula for alpha
    numerator = 4 - n
    denominator = 2 * (n + 8)

    # Calculate alpha
    alpha = (numerator / denominator) * epsilon

    # --- Output the step-by-step calculation ---
    print("This script calculates the specific heat scaling exponent α using the epsilon expansion.")
    print("-" * 70)
    print(f"The calculation is for a system with spatial dimension d = {d}.")
    print(f"We assume the Ising model, which corresponds to n = {n} component(s) of spin.")
    print("\nStep 1: Determine the expansion parameter ϵ (epsilon).")
    print(f"The upper critical dimension for this model is d_c = {d_c}.")
    print(f"ϵ = d_c - d = {d_c} - {d} = {epsilon}")

    print("\nStep 2: Apply the first-order formula for the scaling exponent α.")
    print("The formula is: α = (4 - n) / (2 * (n + 8)) * ϵ")

    print("\nStep 3: Substitute the values and calculate α.")
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = ({numerator}) / (2 * ({n+8})) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    print(f"α = {numerator/6} / {denominator/6}  (simplified fraction)")
    print(f"α = 1 / 6")

    print("\n" + "-" * 70)
    print(f"The final calculated value for the scaling exponent α is: {alpha}")
    print("-" * 70)

if __name__ == '__main__':
    calculate_scaling_exponent_alpha()