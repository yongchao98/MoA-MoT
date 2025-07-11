import math

def calculate_main_diagonal_symmetry():
    """
    Calculates the number of involutions of 8 elements.
    This corresponds to configurations symmetric along the main diagonal.
    The number of involutions a_n follows the recurrence:
    a_n = a_{n-1} + (n-1) * a_{n-2}
    """
    # a[i] will store a_i
    a = [1, 1]  # Base cases: a_0 = 1, a_1 = 1
    for i in range(2, 9):
        a_i = a[i-1] + (i-1) * a[i-2]
        a.append(a_i)
    return a[8]

def calculate_anti_diagonal_symmetry():
    """
    Calculates the number of permutations of 8 elements that commute
    with the permutation s(x) = 9-x.
    This corresponds to configurations symmetric along the anti-diagonal.
    For n=2m, the formula is 2^m * m!.
    """
    m = 4  # For an 8x8 board, n=8, so m=4
    return (2**m) * math.factorial(m)

def calculate_both_diagonals_symmetry():
    """
    Calculates the number of permutations that are both involutions and
    commute with s(x) = 9-x.
    This corresponds to configurations symmetric along both diagonals.
    The number c_{2m} follows the recurrence:
    c_{2m} = 2*c_{2m-2} + 2*(m-1)*c_{2m-4}
    """
    # c[m] will store c_{2m}
    c = [1, 2]  # Base cases: c_0 = 1, c_2 = 2
    for m in range(2, 5):  # We need to calculate up to m=4 (for n=8)
        c_2m = 2 * c[m-1] + 2 * (m-1) * c[m-2]
        c.append(c_2m)
    return c[4]

# Calculate the number of configurations for each case
n_main = calculate_main_diagonal_symmetry()
n_anti = calculate_anti_diagonal_symmetry()
n_both = calculate_both_diagonals_symmetry()

# Apply the Principle of Inclusion-Exclusion
total_configurations = n_main + n_anti - n_both

# Print the final result in the required format
print("The number of configurations symmetric along the main diagonal is: " + str(n_main))
print("The number of configurations symmetric along the anti-diagonal is: " + str(n_anti))
print("The number of configurations symmetric along both diagonals is: " + str(n_both))
print("Using the Principle of Inclusion-Exclusion, the total number of possible configurations is:")
print(f"{n_main} + {n_anti} - {n_both} = {total_configurations}")
<<<1072>>>