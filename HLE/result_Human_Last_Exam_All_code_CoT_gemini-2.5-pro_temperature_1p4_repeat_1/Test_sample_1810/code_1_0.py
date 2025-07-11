import sys

# It is necessary to increase the recursion limit for deep continued fractions
sys.setrecursionlimit(2000)

memo = {}
def N_recursive(coeffs):
    """
    Calculates the numerator of a continued fraction [x_1, x_2, ...]
    using a recursive formula with memoization.
    K_n = a_n * K_{n-1} + K_{n-2}
    """
    coeffs_tuple = tuple(coeffs)
    if coeffs_tuple in memo:
        return memo[coeffs_tuple]
    
    if not coeffs:
        # Numerator of an empty continued fraction K() is defined as 1.
        # This is needed for the k=2 case where N(a_2,...,a_{k-1}) is N().
        return 1
    
    if len(coeffs) == 1:
        return coeffs[0]
    
    # Recursive step based on K_n = a_n * K_{n-1} + K_{n-2}
    # To implement this, we use the equivalent symmetric recurrence:
    # K(x_1,...,x_n) = x_1 * K(x_2,...,x_n) + K(x_3,...,x_n)
    
    # K(x_1...x_n)
    term1 = coeffs[0] * N_recursive(coeffs[1:])
    
    # K(x_3...x_n)
    if len(coeffs) > 2:
        term2 = N_recursive(coeffs[2:])
    else:
        term2 = 1 # K() = 1

    result = term1 + term2
    memo[coeffs_tuple] = result
    return result

def solve_for_ck_demo():
    """
    Demonstrates the identity for a sample k and a_i list, and calculates c_k.
    """
    k = 4
    # Let a_i be a list of positive integers a_1, a_2, ..., a_k
    a = [1, 2, 3, 4] 
    
    print(f"Demonstrating for k = {k} and a = {a}\n")

    # Construct the sequences for the equation
    # LHS: N[a_2,..., a_k+1, a_k,...,a_1]
    # RHS: N[a_1,...,a_k, a_k,...,a_2] + c_k
    
    # Note: Python list indices are 0-based.
    # a_1 is a[0], a_k is a[k-1].
    a_i = [val for val in a]

    # Sequence for the LHS term
    seq_LHS = a_i[1:k-1] + [a_i[k-1] + 1] + a_i[::-1]
    
    # Sequence for the first term on the RHS
    seq_RHS1 = a_i[0:k] + a_i[1:k][::-1]

    # Calculate the values of the numerators
    val_LHS = N_recursive(seq_LHS)
    val_RHS1 = N_recursive(seq_RHS1)

    # Derive c_k using the formula: c_k = N(a_2..a_{k-1}) * N(a_1..a_k)
    # N(a_2..a_{k-1}) -> a_i[1:k-1]
    # N(a_1..a_k) -> a_i[0:k]
    
    term_c1 = N_recursive(a_i[1:k-1])
    term_c2 = N_recursive(a_i[0:k])
    c_k_val = term_c1 * term_c2
    
    # Print the full equation with the calculated values
    print("The equation is:")
    print(f"N{seq_LHS} = N{seq_RHS1} + c_{k}\n")

    print("Numerically, this is:")
    print(f"{val_LHS} = {val_RHS1} + {c_k_val}\n")

    # Verification
    if val_LHS == val_RHS1 + c_k_val:
        print("The identity holds true for this example.")
        print(f"The value of c_{k} is {c_k_val}.")
    else:
        print("There is a discrepancy in the identity.")

if __name__ == "__main__":
    solve_for_ck_demo()
