import math

def riemann_tensor_components_kahler(m):
    """
    Calculates the number of independent real components of the Riemann tensor
    on a Kähler manifold of complex dimension m.

    Args:
        m (int): The complex dimension of the Kähler manifold (must be a positive integer).

    Returns:
        int: The number of independent real components.
    """
    if not isinstance(m, int) or m <= 0:
        raise ValueError("The complex dimension m must be a positive integer.")

    print(f"Calculating for a Kähler manifold of complex dimension m = {m} (real dimension n = {2*m}).")

    if m == 1:
        # For m=1, the manifold is a Riemann surface. Curvature is determined by the scalar curvature.
        # The general formula for a 2-manifold is n^2(n^2-1)/12 = 2^2(3)/12 = 1.
        # Our formula: m^2 = 1^2 = 1.
        num_components = 1
        print("For m = 1, the manifold is a Riemann surface, and the Bochner tensor vanishes.")
        print(f"The number of components is given by m^2.")
        print(f"Result = {m}^2 = {num_components}")

    elif m == 2:
        # For m=2, the manifold is a Kähler surface. The Bochner tensor vanishes identically.
        # The curvature is entirely determined by the Ricci tensor.
        # The number of components in the Ricci tensor (an m x m Hermitian matrix) is m^2.
        num_components = m**2
        print("For m = 2, the Bochner tensor vanishes, which is a special case.")
        print("The curvature is determined by the Ricci tensor, which has m^2 components.")
        print(f"Result = {m}^2 = {num_components}")

    else: # m >= 3
        # For the generic case m >= 3, the formula is m^2(m^2+1)/2.
        # This formula also works for m=1, but not for m=2.
        num_components = (m**2 * (m**2 + 1)) // 2
        print("For m >= 3, the generic formula for the number of components applies.")
        print(f"The formula is m^2 * (m^2 + 1) / 2.")
        print(f"Result = ({m}^2 * ({m}^2 + 1)) / 2 = ({m**2} * {m**2+1}) / 2 = {num_components}")
        
    print("-" * 20)
    return num_components

# Example calculations for a few dimensions.
components_m1 = riemann_tensor_components_kahler(1)
components_m2 = riemann_tensor_components_kahler(2)
components_m3 = riemann_tensor_components_kahler(3)
components_m4 = riemann_tensor_components_kahler(4)

final_answer = f"""
For m=1: {components_m1}
For m=2: {components_m2}
For m=3: {components_m3}
For m=4: {components_m4}
The number is given by a piecewise function: 1 for m=1, m^2 for m=2, and m^2(m^2+1)/2 for m>=3.
"""
# The prompt format is not well-suited for a function as an answer, so I'll output the
# final formula as the answer.
final_formula_str = "The number of components is 1 if m=1, 4 if m=2, and m^2(m^2+1)/2 if m>=3, where m is the complex dimension."
# Since no specific m is given, no single number can be returned.
# The code above provides the user with the tool to find the answer for any m.
