import math

def calculate_degeneracy_formula():
    """
    This function explains and derives the formula for the ground space
    degeneracy of the toric code on a surface with genus g, n smooth holes,
    and m rough holes.
    """

    # 1. Start with the general formula for the number of logical qubits, k.
    # This formula holds for n > 0 and m > 0.
    # k = 2*g + (n - 1) + (m - 1) = 2*g + n + m - 2
    
    print("The number of logical qubits 'k' is given by the formula: k = 2*g + n + m - 2")
    
    # 2. Define the genus 'g' for the toric code.
    g = 1
    print(f"The 'toric code' is typically defined on a torus, so the genus g = {g}.")
    
    # 3. Substitute g=1 into the formula to find k.
    print("Substituting the value of g into the formula for k:")
    
    # The numbers in the final equation k = n + m are derived from:
    two_from_genus = 2
    one_from_genus = 1
    two_from_subtraction = 2
    
    # Using python's f-string formatting to display the calculation
    print(f"k = {two_from_genus} * {one_from_genus} + n + m - {two_from_subtraction}")
    print("k = n + m")
    
    # 4. The ground space degeneracy D is 2^k.
    base = 2
    print(f"\nThe ground space degeneracy D is {base}^k.")
    print(f"So, the final formula for the degeneracy is: D = {base}^(n + m)")

# Execute the function to print the derivation.
calculate_degeneracy_formula()