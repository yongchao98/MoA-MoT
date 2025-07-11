import math

def calculate_min_diameter():
    """
    This function prompts the user for integers n and m, and calculates the
    minimum possible diameter for a tree with n+2 vertices and m leaves.
    """
    try:
        n_str = input("Enter the value for n (a positive integer): ")
        n = int(n_str)
        m_str = input("Enter the value for m (a positive integer): ")
        m = int(m_str)

        if n <= 0 or m <= 0:
            print("\nError: n and m must be positive integers.")
            return

        # A tree with V = n+2 vertices can have at most V-1 = n+1 leaves.
        if m > n + 1:
            print(f"\nError: A tree with n+2 = {n+2} vertices cannot have m = {m} leaves.")
            print(f"The number of leaves (m) must be less than or equal to n+1 = {n+1}.")
            return

        # The minimum diameter is calculated using the formula: ceil(2 * (n + 1) / m)
        numerator = 2 * (n + 1)
        division_result = numerator / m
        min_diameter = math.ceil(division_result)

        print("\n--- Calculation ---")
        print(f"For a tree with {n+2} vertices and {m} leaves, the minimum possible diameter is found using the formula:")
        print("D_min = ceil(2 * (n + 1) / m)")
        print("\nSubstituting the given values:")
        print(f"D_min = ceil(2 * ({n} + 1) / {m})")
        print(f"      = ceil(2 * {n + 1} / {m})")
        print(f"      = ceil({numerator} / {m})")
        print(f"      = ceil({division_result:.4f})")
        print(f"      = {min_diameter}")

    except ValueError:
        print("\nError: Invalid input. Please enter valid integers for n and m.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    calculate_min_diameter()