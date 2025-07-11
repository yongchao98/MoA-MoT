import sys

def find_ideal_ratio():
    """
    This function presents the ideal Ni/Ce ratio for catalytic applications
    based on findings from scientific literature.
    """
    # According to extensive research, a relatively low Ni content is optimal.
    # This maximizes the active interface between Ni and Ceria (the Ni-O-Ce sites)
    # and prevents the sintering (agglomeration) of Ni particles at high temperatures,
    # which would deactivate the catalyst.
    # The ideal ratio is often expressed as a molar ratio.

    ni_molar_part = 1
    ce_molar_part = 9

    # The task requires printing each number in the final "equation" or result.
    print("--- Determining Optimal Ni/Ce Molar Ratio ---")
    print(f"Component 1: Nickel (Ni) molar part = {ni_molar_part}")
    print(f"Component 2: Cerium (Ce) molar part = {ce_molar_part}")
    print("-" * 45)

    # Print the final conclusion based on these values.
    print("Conclusion:")
    print(f"The ideal molar ratio of Ni to Ce to maximize catalytic performance for")
    print(f"the Water Gas Shift and Water Splitting reactions is approximately {ni_molar_part}:{ce_molar_part}.")
    print("\nThis ratio (corresponding to ~10 at% Ni) creates a high population of active sites at the Ni-Ceria interface while maintaining excellent thermal stability.")

# Execute the function to find and print the ratio.
find_ideal_ratio()
