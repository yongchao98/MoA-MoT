import numpy as np

def calculate_l_of_n(n):
    """
    Calculates the value of l(n) based on the derived symbolic formula.

    The symbolic derivation of l(n) = f^(1)(f_M^(2)(f^(3)(P))) results in the formula:
    l(n) = (2*n**2 + 2 - (4*n - 2) * sqrt(n**2 - n + 1)) / n**2

    This function calculates the result for a given n and prints the components
    of the formula as requested.
    """
    if not isinstance(n, int) or n < 5:
        print("Error: n must be an integer greater than or equal to 5.")
        return

    print(f"Calculating l(n) for n = {n}")
    print("The final formula is: l(n) = (A - B * sqrt(C)) / D")
    print("where:")

    # Calculate and print each component of the formula
    A = 2 * n**2 + 2
    print(f"A = 2*n^2 + 2 = 2*{n**2} + 2 = {A}")

    B = 4 * n - 2
    print(f"B = 4*n - 2 = 4*{n} - 2 = {B}")

    C = n**2 - n + 1
    print(f"C = n^2 - n + 1 = {n**2} - {n} + 1 = {C}")

    D = n**2
    print(f"D = n^2 = {n**2} = {D}")

    # Calculate the final result
    result = (A - B * np.sqrt(C)) / D

    # Print the final equation with the computed numbers
    print("\nSubstituting these values into the formula:")
    print(f"l({n}) = ({A} - {B} * sqrt({C})) / {D}")
    print(f"\nThe exact numerical value is: {result}")


if __name__ == '__main__':
    # The problem is defined for n >= 5. Let's use n=5 as an example.
    n_value = 5
    calculate_l_of_n(n_value)