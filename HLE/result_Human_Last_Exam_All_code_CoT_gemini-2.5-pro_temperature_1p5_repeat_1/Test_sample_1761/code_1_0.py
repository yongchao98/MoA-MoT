import math

def h0_Ok(n, k):
    """
    Calculates the dimension of the space of global sections of O(k) on P^n.
    h^0(P^n, O(k)) = C(n+k, k) for k >= 0.
    """
    if k < 0:
        return 0
    # C(n, k) = n! / (k! * (n-k)!)
    # C(n+k, k) = (n+k)! / (k! * n!)
    if k == 0:
        return 1
    if n == 0 and k > 0:
        return 1 # P^0 is a point, O(k) is trivial bundle
    return math.comb(n + k, k)

def complex_dimension_of_global_sections(n):
    """
    Calculates the complex dimension of the space of global sections for
    Omega^1_{P^n} tensored with O(2).
    The result is h^0(P^n, Omega^1(2)).
    """
    if n < 0:
      raise ValueError("n must be a non-negative integer.")
      
    print(f"For n = {n}:")

    # Step 1: Calculate the dimension of the middle terms in the exact sequence.
    # h^0(P^n, O(1))^(n+1)
    h0_O1 = h0_Ok(n, 1)
    dim_middle_term1 = (n + 1) * h0_O1

    # h^0(P^n, O(2))
    dim_middle_term2 = h0_Ok(n, 2)
    
    print("The dimension is derived from the long exact sequence in cohomology associated with:")
    print("0 -> Omega^1(2) -> O(1)^(n+1) -> O(2) -> 0")
    print("\nThis yields the equation for the dimension h^0(Omega^1(2)):")
    print("h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))")
    
    print("\nCalculating the terms:")
    print(f"h^0(P^{n}, O(1)) = C({n}+1, 1) = {h0_O1}")
    print(f"h^0(P^{n}, O(1))^({n}+1) = {n+1} * {h0_O1} = {dim_middle_term1}")
    print(f"h^0(P^{n}, O(2)) = C({n}+2, 2) = {dim_middle_term2}")

    # Step 2: Calculate the final dimension.
    result = dim_middle_term1 - dim_middle_term2
    
    print("\nThe final dimension is:")
    print(f"h^0(Omega^1(2)) = {dim_middle_term1} - {dim_middle_term2} = {result}")

    # For verification, we can also directly use the formula n(n+1)/2
    verification = n * (n + 1) // 2
    print(f"\nVerification with the formula n(n+1)/2 = {n}({n}+1)/2 = {verification}")

    return result

if __name__ == '__main__':
    # You can change the value of n here.
    # The question is for a general n, here we demonstrate with n=3.
    n_value = 3
    final_dimension = complex_dimension_of_global_sections(n_value)
    # The final answer format is specified to be just the content.
    # Since the question is for general n, the content is the formula.
    # For a specific n, e.g. n=3, the answer is a number. Let's return the number for n=3.
    # print(f"<<<{final_dimension}>>>")