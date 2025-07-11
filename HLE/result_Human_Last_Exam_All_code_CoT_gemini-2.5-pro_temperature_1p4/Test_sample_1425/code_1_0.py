import math

def display_partition_function():
    """
    This function prints the symbolic formula for the partition function Z
    for a system with Hamiltonian H = -μN in the grand canonical ensemble.
    The formula is derived assuming the particles are bosons.
    """

    # Using unicode characters for better readability
    beta = "\u03B2"  # β
    mu = "\u03BC"    # μ
    prod_symbol = "\u220F"  # Π
    sum_symbol = "\u2211"   # Σ

    print("The grand canonical partition function Z for a system of bosons with Hamiltonian H = -μN is:")
    # The final equation includes the numbers 1 and 2 explicitly as requested.
    print(f"Z = {prod_symbol}_k [ 1 / (1 - exp(2 * {beta} * {mu})) ]")
    print("\nThis result holds for μ < 0 to ensure the series converges.")
    print("\nThe logarithm of the partition function, which is often used to find thermodynamic potentials, is:")
    print(f"log(Z) = -{sum_symbol}_k log( 1 - exp(2 * {beta} * {mu}) )")

if __name__ == '__main__':
    display_partition_function()