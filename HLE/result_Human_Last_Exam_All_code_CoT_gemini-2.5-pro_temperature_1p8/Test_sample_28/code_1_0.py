import sys

def compute_vc_dimension(z, T):
    """
    Computes the VC dimension for the hypothesis class H_{z-ones}.

    H_{z-ones} is the set of all hypotheses h: X -> {0,1} on a domain X of size T
    such that the number of points x for which h(x)=1 is exactly z.

    Args:
        z (int): The required number of '1' labels. Must be a non-negative integer.
        T (int): The total size of the domain X. Must be a non-negative integer.
    """
    print(f"----- Calculating VC Dimension for z={z} and T={T} -----")

    if not isinstance(z, int) or not isinstance(T, int) or T < 0 or z < 0:
        print("Error: z and T must be non-negative integers.")
        return

    # Case 1: The hypothesis class is empty because it's impossible to choose
    # z elements from a set of size T.
    if z > T:
        print("The hypothesis class is empty because z > T.")
        print("VC Dimension = -1")
        return

    # Case 2: The hypothesis class is non-empty. The VC dimension is min(z, T-z).
    result = min(z, T - z)
    # The user requested printing the full equation
    print(f"VC Dimension = min(z, T - z)")
    print(f"             = min({z}, {T} - {z})")
    print(f"             = min({z}, {T - z})")
    print(f"             = {result}")

if __name__ == '__main__':
    # Example 1: A standard case
    compute_vc_dimension(z=10, T=50)

    # Example 2: An edge case where z is close to T
    compute_vc_dimension(z=48, T=50)
    
    # Example 3: An edge case where z equals T
    compute_vc_dimension(z=50, T=50)
    
    # Example 4: A case where the hypothesis class is empty
    compute_vc_dimension(z=15, T=10)

    # Example 5: Another standard case
    compute_vc_dimension(z=3, T=10)