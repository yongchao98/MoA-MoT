import numpy as np
from itertools import permutations

def get_permutation_class(p_tuple):
    """
    Assigns a permutation of S3 to its conjugacy class.
    Class 0: Identity
    Class 1: Transpositions (3 elements)
    Class 2: 3-cycles (2 elements)
    """
    if p_tuple == (0, 1, 2):
        return 0  # Identity
    
    fixed_points = sum(1 for i in range(3) if p_tuple[i] == i)
    if fixed_points == 1:
        return 1  # Transposition
    if fixed_points == 0:
        return 2  # 3-cycle
    # This part should not be reached for S3 permutations
    return -1

def calculate_immanant_verbose(M, chi_name, chi_values):
    """
    Calculates the immanant for a given matrix M and character chi,
    printing the full calculation.
    """
    n = M.shape[0]
    if n != 3:
        raise ValueError("This implementation is specific to n=3.")

    print(f"\nCalculating immanant for character: {chi_name}")
    print("Formula: Immanant = sum_{p in S3} [chi(p) * product(M[i, p(i)])]")
    
    total_immanant = 0
    calculation_steps = []
    
    # Iterate over all permutations in S3
    for p in permutations(range(n)):
        term_product = 1
        product_factors = []
        for i in range(n):
            val = M[i, p[i]]
            term_product *= val
            product_factors.append(str(val))

        # Get character value based on permutation class
        perm_class = get_permutation_class(p)
        char_val = chi_values[perm_class]
        
        term = char_val * term_product
        total_immanant += term
        
        # Format the string for this step of the calculation
        # e.g., "(2) * ((-5) * (3) * (2))"
        step_str = f"({char_val}) * ({' * '.join(product_factors)})"
        calculation_steps.append(step_str)

    print("Calculation: " + " + ".join(calculation_steps))
    print(f"Result: {total_immanant}")
    return total_immanant

# For n=3, a "Mercer matrix" M_3 that satisfies the problem's conditions is:
M3 = np.array([
    [-5, 1, -10],
    [-1, 3, -2],
    [1, 4, 2]
])

print("For the case n=3, the specific matrix M_3 chosen is:")
print(M3)

# The irreducible characters of S3 are defined by their values on the
# conjugacy classes (identity, transposition, 3-cycle).
characters_s3 = {
    "Trivial (Permanent)": [1, 1, 1],
    "Sign (Determinant)": [1, -1, 1],
    "Standard": [2, 0, -1]
}

immanant_values = {}
# Calculate the immanant for each character
for name, values in characters_s3.items():
    immanant_values[name] = calculate_immanant_verbose(M3, name, values)

# Find the largest immanant in absolute value
largest_immanant = 0
for val in immanant_values.values():
    if abs(val) > largest_immanant:
        largest_immanant = abs(val)

print("\n--- Summary ---")
print(f"The calculated immanants are: {immanant_values}")
print(f"The largest immanant in absolute value is: {largest_immanant}")
