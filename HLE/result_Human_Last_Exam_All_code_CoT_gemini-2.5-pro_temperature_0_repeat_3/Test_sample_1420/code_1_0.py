def calculate_alpha_exponent():
    """
    Calculates the specific heat scaling exponent alpha using the first-order
    epsilon expansion for the 3D Ising model.
    """
    # Step 1 & 2: Define the dimensions and model parameter n
    d_c = 4  # Upper critical dimension for the O(n) model
    d = 3    # Spatial dimension
    n = 1    # Number of components for the order parameter (n=1 for Ising model)

    # Step 3: Calculate epsilon
    epsilon = d_c - d

    # Step 4 & 5: Use the first-order epsilon expansion formula for alpha
    # Formula: alpha = (4 - n) / (n + 8) * epsilon
    numerator = 4 - n
    denominator = n + 8

    # Step 6: Calculate the final value of alpha
    alpha = (numerator / denominator) * epsilon

    # Print the explanation and the calculation step-by-step
    print("Calculation of the specific heat scaling exponent α using the first-order epsilon expansion:")
    print(f"The upper critical dimension is d_c = {d_c}.")
    print(f"The spatial dimension is d = {d}.")
    print(f"The expansion parameter is ϵ = d_c - d = {d_c} - {d} = {epsilon}.")
    print(f"For the Ising universality class, the number of order parameter components is n = {n}.")
    print("\nThe first-order expansion formula for α is: α = (4 - n) / (n + 8) * ϵ")
    print("\nSubstituting the values into the equation:")
    print(f"α = ({4} - {n}) / ({n} + {8}) * {epsilon}")
    print(f"α = {numerator} / {denominator} * {epsilon}")
    print(f"α = {alpha}")

calculate_alpha_exponent()