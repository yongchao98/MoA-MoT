import math

def solve_csp_counting():
    """
    Calculates the number of NP-hard sets S for the given constraint satisfaction problem.
    """
    # 1. Total number of subsets S
    base_set_size = 5
    total_sets = 2**base_set_size

    print(f"Let S be a subset of {{0, 1, 2, 3, 4}}.")
    print(f"The total number of possible sets S is 2^{base_set_size} = {total_sets}.")
    print("-" * 30)

    # 2. Count the P-time cases
    print("A problem is in P if its constraint set S is 0-valid, 1-valid, Horn, dual-Horn, or affine.")
    print("We will count the number of such sets S.\n")

    # Count sets that are 0-valid or 1-valid
    # Number of sets containing a specific element (e.g., 0) is 2^(n-1)
    num_0_valid = 2**(base_set_size - 1)
    # Number of sets containing another specific element (e.g., 4) is 2^(n-1)
    num_4_valid = 2**(base_set_size - 1)
    # Number of sets containing both 0 and 4 is 2^(n-2)
    num_0_and_4_valid = 2**(base_set_size - 2)
    
    # Using the Principle of Inclusion-Exclusion for |A U B| = |A| + |B| - |A intersect B|
    num_0_or_4_valid = num_0_valid + num_4_valid - num_0_and_4_valid

    print("Counting P-time sets:")
    print(f"Number of 0-valid sets (S contains 0): 2^({base_set_size - 1}) = {num_0_valid}")
    print(f"Number of 1-valid sets (S contains 4): 2^({base_set_size - 1}) = {num_4_valid}")
    print(f"Number of sets containing both 0 and 4: 2^({base_set_size - 2}) = {num_0_and_4_valid}")
    print(f"Number of sets that are 0-valid OR 1-valid = {num_0_valid} + {num_4_valid} - {num_0_and_4_valid} = {num_0_or_4_valid}\n")

    # Consider other P-time cases not covered above
    # These are sets S that do NOT contain 0 and do NOT contain 4.
    # Such sets are subsets of {1, 2, 3}.
    # We check which of these are Horn, dual-Horn, or affine.
    # - Horn sets not containing 0 or 4: S = {}
    # - Dual-Horn sets not containing 0 or 4: S = {}
    # - Affine sets not containing 0 or 4: S = {} and S = {1, 3}
    num_additional_p_cases = 2 # For S={} and S={1,3}
    print(f"We check for remaining P-time cases among subsets of {{1, 2, 3}}.")
    print(f"These are the cases that are neither 0-valid nor 1-valid.")
    print(f"The additional P-time sets found are S={{}} and S={{1, 3}}.")
    print(f"Number of additional P-time sets = {num_additional_p_cases}\n")

    # Total P-time sets
    total_p_time = num_0_or_4_valid + num_additional_p_cases
    print(f"Total number of P-time sets = {num_0_or_4_valid} + {num_additional_p_cases} = {total_p_time}")
    print("-" * 30)

    # 3. Calculate NP-hard cases
    num_np_hard = total_sets - total_p_time
    print("The number of NP-hard sets is the total number of sets minus the P-time sets.")
    print(f"Number of NP-hard sets = {total_sets} - {total_p_time} = {num_np_hard}")

if __name__ == '__main__':
    solve_csp_counting()