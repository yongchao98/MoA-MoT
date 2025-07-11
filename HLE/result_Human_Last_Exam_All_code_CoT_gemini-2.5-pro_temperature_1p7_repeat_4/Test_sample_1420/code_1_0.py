def calculate_specific_heat_exponent():
    """
    Calculates the specific heat exponent alpha (α) using the first-order
    epsilon expansion for a system in d=3.
    """
    # --- Step 1: Define parameters ---
    # Spatial dimension
    d = 3
    # Upper critical dimension for the phi^4 model
    d_c = 4
    # Number of order parameter components (n=1 for the Ising model)
    n = 1

    # --- Step 2: Calculate epsilon ---
    epsilon = d_c - d

    print(f"Calculating the specific heat exponent α for a system with spatial dimension d = {d}.")
    print(f"The Ising model universality class (n={n}) is assumed.")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The expansion parameter is ε = d_c - d = {d_c} - {d} = {epsilon}.")
    print("-" * 30)

    # --- Step 3: Apply the formula ---
    print("The first-order formula for α is: (4 - n) / (2 * (n + 8)) * ε")
    print("Substituting the values into the equation:")

    numerator = 4 - n
    denominator_part1 = n + 8
    denominator = 2 * denominator_part1
    
    # Calculate the final value of alpha
    alpha = numerator / denominator * epsilon

    # --- Step 4: Print the calculation steps and result ---
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = {numerator} / (2 * {denominator_part1})")
    print(f"α = {numerator} / {denominator}")
    print(f"α = {alpha}")

if __name__ == "__main__":
    calculate_specific_heat_exponent()
