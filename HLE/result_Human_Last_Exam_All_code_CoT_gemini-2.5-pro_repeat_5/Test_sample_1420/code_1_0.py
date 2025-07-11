def calculate_alpha_epsilon_expansion():
    """
    Calculates the specific heat scaling exponent alpha using the first-order
    epsilon expansion for the Ising model in d=3.
    """
    # 1. Define the parameters
    d = 3  # Spatial dimension
    d_c = 4  # Upper critical dimension for the O(n) model
    n = 1  # Number of components of the order parameter (n=1 for Ising model)

    # 2. Calculate epsilon
    epsilon = d_c - d

    # 3. Perform the calculation using the first-order formula for alpha
    # Formula: alpha = (4 - n) / (2 * (n + 8)) * epsilon
    numerator = 4 - n
    denominator = 2 * (n + 8)
    alpha = (numerator / denominator) * epsilon

    # 4. Print the explanation and result
    print("Calculation of the specific heat scaling exponent α using the first-order epsilon expansion.")
    print("-" * 80)
    print(f"The calculation is performed for a system with spatial dimension d = {d}.")
    print(f"The relevant upper critical dimension (d_c) is {d_c}.")
    print(f"We assume the Ising model, where the order parameter has n = {n} component.")
    print("\nStep 1: Calculate the expansion parameter epsilon (ϵ)")
    print(f"ϵ = d_c - d")
    print(f"ϵ = {d_c} - {d} = {epsilon}")
    print("\nStep 2: Use the first-order expansion formula for α")
    print("α = (4 - n) / (2 * (n + 8)) * ϵ")
    print("\nStep 3: Substitute the values into the formula")
    print(f"α = (4 - {n}) / (2 * ({n} + 8)) * {epsilon}")
    print(f"α = {numerator} / (2 * {n + 8}) * {epsilon}")
    print(f"α = {numerator} / {2 * (n + 8)}")
    print(f"α = {alpha}")

# Execute the function to print the output
calculate_alpha_epsilon_expansion()