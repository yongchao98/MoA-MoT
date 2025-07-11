import math

def calculate_toric_code_gsd(n, m):
    """
    Calculates the ground space degeneracy (GSD) of the toric code
    on a sphere with n smooth holes and m rough holes.

    The calculation is based on the formula for the case where both
    n > 0 and m > 0.

    Args:
        n (int): The number of smooth holes.
        m (int): The number of rough holes.

    Returns:
        A tuple containing the number of logical qubits and the GSD,
        or an error message string if n or m are not positive.
    """
    if n <= 0 or m <= 0:
        return (None, "This formula assumes n > 0 and m > 0. The physics changes for cases where one type of hole is absent.")

    # Number of logical qubits k = n + m - 2
    num_qubits = n + m - 2

    # Ground Space Degeneracy = 2^k
    gsd = int(math.pow(2, num_qubits))

    return (num_qubits, gsd)

def main():
    """
    Main function to demonstrate the calculation with example values.
    """
    print("This program calculates the ground space degeneracy (GSD) of the toric code")
    print("on a surface with n smooth holes and m rough holes.")
    print("The governing formula for cases where n > 0 and m > 0 is: GSD = 2^(n + m - 2)")
    print("-" * 20)

    # Example values for n and m
    n_example = 3
    m_example = 4

    print(f"Let's calculate the GSD for an example with n = {n_example} and m = {m_example}.\n")

    num_qubits, gsd = calculate_toric_code_gsd(n_example, m_example)

    if gsd is not None:
        print("The final equation is derived as follows:")
        # The prompt requires outputting each number in the final equation.
        print(f"GSD = 2^({n_example} + {m_example} - 2)")
        print(f"GSD = 2^({num_qubits})")
        print(f"GSD = {gsd}")

if __name__ == "__main__":
    main()
