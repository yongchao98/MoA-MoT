import math

def calculate_and_print_destabilizers(n):
    """
    Calculates and prints the number of different destabilizer sets for a given number of qubits n.

    The stabilizer set is {Z_1, ..., Z_n}.
    The formula for the total number of destabilizer sets is:
    N = (phase_choices) * (XY_choices) * (Z_choices)
      = 4^n * 2^n * 2^(n*(n-1)/2)
      = 2^((n^2 + 5n)/2)

    Args:
        n (int): The number of qubits, must be a non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Error: The number of qubits 'n' must be a non-negative integer.")
        return

    # 1. Calculate the number of choices for phases (4 options: +1, -1, +i, -i for each of n destabilizers)
    # Using 4**n might lead to very large numbers, so we work with powers of 2.
    # 4^n = (2^2)^n = 2^(2*n)
    power_from_phases = 2 * n

    # 2. Calculate the number of choices for the anti-commuting part (X or Y for each of n destabilizers)
    # This gives 2^n choices.
    power_from_xy = n

    # 3. Calculate the number of choices for the commuting Z-parts
    # This is determined by a symmetric n x n matrix with zero diagonal, which has n*(n-1)/2 independent entries.
    # This gives 2^(n*(n-1)/2) choices.
    # Using integer division // ensures the result is an integer.
    power_from_z = n * (n - 1) // 2
    
    # 4. Total number is the product of these choices. In terms of powers of 2, we sum the exponents.
    total_power = power_from_phases + power_from_xy + power_from_z
    
    # Calculate the components of the final equation for printing
    num_phases = 4**n
    num_xy = 2**n
    num_z = 2**power_from_z
    total_sets = num_phases * num_xy * num_z

    print(f"For n = {n} qubits:")
    print("The total number of different destabilizer sets is calculated as:")
    print("(Number of phase choices) * (Number of X/Y choices) * (Number of Z-part choices)")
    
    # Print the final equation with the calculated numbers
    print(f"\nFinal Equation: {num_phases} * {num_xy} * {num_z} = {total_sets}")
    
    # Optional: Print the result as a power of 2 for verification
    if total_power == (n**2 + 5*n)/2:
        print(f"This is equal to 2^{total_power}.")
    print("-" * 30)

if __name__ == '__main__':
    # Example calculations for different values of n
    calculate_and_print_destabilizers(1)
    calculate_and_print_destabilizers(2)
    calculate_and_print_destabilizers(3)
    calculate_and_print_destabilizers(5)
