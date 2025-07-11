import sys

def get_pi_so(i, m):
    """
    Calculates the i-th homotopy group of SO(m).
    This function provides values in the stable range and for pi_1.
    """
    if i <= 0:
        return "0"
    
    # pi_1(SO(m))
    if i == 1:
        if m == 1:
            return "0" # SO(1) is a point
        if m == 2:
            return "Z" # Z for integers
        if m >= 3:
            return "Z_2" # Z_2 for integers mod 2

    # Stable range check for i >= 2
    # The stable range is roughly i < m-1.
    if i >= m - 1:
        return f"pi_{i}(SO({m})) (unstable, complex)"

    # Bott periodicity for stable homotopy groups of SO for i >= 2
    # pi_i(SO) is periodic with period 8.
    # We use the isomorphism pi_i(SO) = pi_{i+1}(O) for i>=1
    # and the stable groups of O are pi_k(O) with k=i+1
    k = i
    rem = k % 8
    if rem == 3 or rem == 7:
        return "Z"
    elif rem == 0 or rem == 1:
        return "Z_2"
    else: # rem in {2, 4, 5, 6}
        return "0"

def main():
    """
    Main function to identify and calculate the obstruction groups.
    """
    # You can change these values
    n = 4  # Dimension of the sphere is n-1
    k = 5  # Rank of the vector bundle is 2k

    if n < 1 or k < 1:
        print("Please provide positive integers for n and k.")
        return

    m = 2 * k

    print("The homotopy-theoretic obstructions for the two paths to be homotopic are classified by elements of two groups.")
    print("-" * 30)

    # First obstruction group
    group1_name = f"pi_1(SO({m}))"
    group1_value = get_pi_so(1, m)
    print(f"1. The first obstruction group is {group1_name}.")
    print(f"   For the given k={k} (m={m}), this group is: {group1_value}")
    
    print("-" * 30)

    # Second obstruction group
    i = n + 1
    group2_name = f"pi_{i}(SO({m}))"
    group2_value = get_pi_so(i, m)
    print(f"2. The second obstruction group is {group2_name}.")
    print(f"   For the given n={n}, k={k} (i={i}, m={m}), this group is: {group2_value}")
    print("-" * 30)
    
    print("\nNote:")
    print("Z represents the group of integers.")
    print("Z_2 represents the cyclic group of order 2.")
    print("0 represents the trivial group.")
    print("The calculation for higher homotopy groups is provided within the 'stable range'.")
    print("Outside this range, the groups are 'unstable' and generally more complex to compute.")


if __name__ == "__main__":
    main()
