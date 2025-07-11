def calculate_riemann_components_kahler(complex_dimension_n):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of a given complex dimension n.

    The function also prints the step-by-step derivation of the result.
    """
    if not isinstance(complex_dimension_n, int) or complex_dimension_n <= 0:
        print("Error: Please provide a positive integer for the complex dimension.")
        return

    print(f"Calculating for a Kähler manifold of complex dimension n = {complex_dimension_n}.")
    print("-----------------------------------------------------------------")
    print("Step 1: The tensor is described by an N x N Hermitian matrix.")
    print("        N is the dimension of the space of symmetric 2-tensors.")
    
    # Calculate N
    n = complex_dimension_n
    n_plus_1 = n + 1
    numerator_N = n * n_plus_1
    N = numerator_N // 2

    print(f"\nStep 2: Calculate N for n = {n}.")
    print(f"        N = (n * (n + 1)) / 2")
    print(f"        N = ({n} * {n_plus_1}) / 2 = {numerator_N} / 2 = {N}")

    # Calculate the final result
    result = N ** 2

    print("\nStep 3: Calculate the total number of independent components.")
    print(f"        The number of independent real components is N^2.")
    print(f"        Components = {N}^2 = {result}")

    print("\n--- Final Equation Breakdown ---")
    print("The final formula for the number of components is ((n * (n + 1)) / 2)^2.")
    print(f"For n = {n}:")
    # This print statement shows each number in the final equation as requested.
    print(f"Components = (({n} * ({n} + 1)) / 2)^2 = (({n} * {n_plus_1}) / 2)^2 = ({numerator_N} / 2)^2 = ({N})^2 = {result}")
    print("-----------------------------------------------------------------")

# You can change the complex dimension n here.
# For example, for a K3 surface, n = 2. For a Calabi-Yau threefold, n = 3.
# Let's use n=2 as an example.
complex_dimension = 2
calculate_riemann_components_kahler(complex_dimension)

# As another example, let's try n=3
complex_dimension_three = 3
calculate_riemann_components_kahler(complex_dimension_three)