import math

def calculate_lower_bound(q):
    """
    Calculates the exact lower bound m >= (1/(q+1)) * C(2q, q).
    
    Args:
        q (int): The sparsity parameter.
    """
    if not isinstance(q, int) or q < 1:
        print("Error: q must be a positive integer.")
        return

    # Calculate combination C(2q, q)
    # C(n, k) = n! / (k! * (n-k)!)
    # C(2q, q) = (2q)! / (q! * q!)
    try:
        comb = math.comb(2 * q, q)
    except ValueError:
        print(f"Error: Could not compute combinations for q={q}. It might be too large.")
        return
        
    lower_bound = comb / (q + 1)
    
    # The question asks for the formula. We print the components of the formula and the result.
    # We output the components of the exact formula derived.
    # The equation is m >= C(2q,q) - C(2q, q-1)
    
    comb_minus_1 = math.comb(2 * q, q - 1)
    
    print(f"The lower bound for m is derived from the formula: C(2q, q) - C(2q, q-1)")
    print(f"For q = {q}:")
    print(f"C({2*q}, {q}) = {comb}")
    print(f"C({2*q}, {q-1}) = {comb_minus_1}")
    print(f"The calculated lower bound is {comb} - {comb_minus_1} = {lower_bound}")
    print(f"This can also be written as (1/({q}+1)) * C({2*q}, {q})")
    print(f"Which is (1/{q+1}) * {comb} = {lower_bound}")


# Example usage with a value for q.
# You can change this value to see the bound for different q.
q_value = 4
calculate_lower_bound(q_value)