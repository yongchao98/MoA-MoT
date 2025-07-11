import sys

def calculate_riemann_components_kahler():
    """
    Calculates and explains the number of independent components of the Riemann
    tensor on a Kähler manifold of a given complex dimension.
    """
    try:
        # Prompt the user for the complex dimension n
        n_str = input("Enter the complex dimension 'n' of the Kähler manifold: ")
        n = int(n_str)
        if n < 0:
            print("Dimension cannot be negative.", file=sys.stderr)
            return
    except ValueError:
        print("Invalid input. Please enter a non-negative integer.", file=sys.stderr)
        return

    # Step 1: Calculate N, the size of the equivalent Hermitian matrix.
    # N is the number of symmetric pairs from n elements.
    # The formula is N = n * (n + 1) / 2
    n_plus_1 = n + 1
    N = n * n_plus_1 / 2
    # Ensure N is an integer for clarity in printing, as it always will be.
    N_int = int(N)

    # Step 2: Calculate the number of independent components, which is N^2.
    result = N_int ** 2

    # Print the explanation of the calculation with the actual numbers.
    print(f"\nFor a Kähler manifold of complex dimension n = {n}:")
    print("\nThe formula for the number of independent components is (n * (n + 1) / 2)^2.")
    
    print("\nStep 1: Calculate the intermediate value N = n * (n + 1) / 2")
    print(f"N = {n} * ({n} + 1) / 2")
    print(f"N = {n} * {n_plus_1} / 2 = {N_int}")

    print("\nStep 2: The number of independent components is N^2.")
    print(f"Number of components = {N_int}^2 = {result}")

    print(f"\nFinal Answer: The Riemann tensor on a Kähler manifold of complex dimension {n} (real dimension {2*n}) has {result} independent entries.")

if __name__ == '__main__':
    calculate_riemann_components_kahler()