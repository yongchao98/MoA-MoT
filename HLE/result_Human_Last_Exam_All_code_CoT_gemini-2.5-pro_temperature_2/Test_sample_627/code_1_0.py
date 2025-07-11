import re

def get_z_degree(poly_str):
    """
    Finds the highest power of 'z' in a polynomial string.
    Assumes format like 'z**4 - 2*z**2' or a constant like '-1'.
    """
    if 'z' not in poly_str:
        return 0
    # Find all occurrences of z**<number> or z
    powers = re.findall(r'z(?:\*\*(\d+))?', poly_str)
    # Handle standalone 'z' which has a power of 1
    # Convert found powers (as strings) to integers. If a power is not specified (just 'z'), default to 1.
    int_powers = [int(p) if p else 1 for p in powers]
    return max(int_powers)

# HOMFLY-PT polynomial for the three-twist knot (6_1) expressed as a dictionary
# where keys are powers of 'a' and values are coefficient polynomials in 'z'.
# P(a, z) = (-1)a^6 + (-2 + z^2)a^4 + (z^4 - 2*z^2)a^2
homfly_poly_in_a = {
    6: "-1",
    4: "-2 + z**2",
    2: "z**4 - 2*z**2"
}

# Identify min and max powers of 'a'
a_max = max(homfly_poly_in_a.keys())
a_min = min(homfly_poly_in_a.keys())

print(f"The HOMFLY-PT polynomial P(a,z) is rearranged by powers of 'a':")
print(f"P(a,z) = ({homfly_poly_in_a[6]})*a**6 + ({homfly_poly_in_a[4]})*a**4 + ({homfly_poly_in_a[2]})*a**2\n")

# Get the corresponding coefficient polynomials
P_a_max_str = homfly_poly_in_a[a_max]
P_a_min_str = homfly_poly_in_a[a_min]

# Calculate the degree in z for each coefficient polynomial
deg_z_P_a_max = get_z_degree(P_a_max_str)
deg_z_P_a_min = get_z_degree(P_a_min_str)

print(f"The minimum power of 'a' is a_min = {a_min}.")
print(f"The coefficient of a**{a_min} is P_amin(z) = {P_a_min_str}.")
print(f"The degree of P_amin(z) is deg_z(P_amin) = {deg_z_P_a_min}.\n")

print(f"The maximum power of 'a' is a_max = {a_max}.")
print(f"The coefficient of a**{a_max} is P_amax(z) = {P_a_max_str}.")
print(f"The degree of P_amax(z) is deg_z(P_amax) = {deg_z_P_a_max}.\n")

# Apply Vogel's algorithm formula
# b(K) <= 1 + (1/2) * (deg_z(P_amin) - deg_z(P_amax))
upper_bound = 1 + 0.5 * (deg_z_P_a_min - deg_z_P_a_max)

print("Using Vogel's algorithm formula: b(K) <= 1 + (1/2) * (deg_z(P_amin) - deg_z(P_amax))")
# Final print statement as requested, showing each number in the equation.
print(f"Upper Bound = 1 + (1/2) * ({deg_z_P_a_min} - {deg_z_P_a_max}) = {int(upper_bound)}")
