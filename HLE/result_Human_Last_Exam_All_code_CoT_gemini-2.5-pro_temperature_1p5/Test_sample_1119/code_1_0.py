import math

def calculate_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on a combinatorial formula
    derived using the Principle of Inclusion-Exclusion.

    The final formula is:
    Count = Sum_{p=0 to K-1} [(-1)^p * C(K-1, p) * C(N - p*M, K)]
    
    where C(n, k) is the binomial coefficient "n choose k".
    
    This function calculates this sum and also formats a string to show the
    step-by-step calculation.
    """
    
    total = 0
    symbolic_parts = []
    value_parts = []
    
    # The sum goes from p = 0 up to K-1
    for p in range(K):
        
        # The term involves C(N - p*M, K). If N - p*M < K, this is 0.
        # All subsequent terms will also be 0, so we can stop early.
        if N - p * M < K:
            break
            
        # Calculate the binomial coefficients for the current term
        comb1 = math.comb(K - 1, p)
        comb2 = math.comb(N - p * M, K)
        
        # The value of the p-th term in the sum (without the sign)
        term_value = comb1 * comb2
        
        # Determine the sign of the term based on p (from (-1)^p)
        if p % 2 == 0:
            # For even p, add to the total. The sign is '+' (or none for the first term).
            total += term_value
            sign_str = " + " if p > 0 else ""
        else:
            # For odd p, subtract from the total. The sign is '-'.
            total -= term_value
            sign_str = " - "

        # Build the symbolic and numeric parts of the equation string for this term
        symbolic_parts.append(f"{sign_str}C({K-1}, {p})*C({N - p*M}, {K})")
        value_parts.append(f"{sign_str}{term_value}")
        
    # Assemble the final equation string for printing
    # Example format: C(2,0)*C(10,3) - C(2,1)*C(6,3) = 120 - 40 = 80
    equation_str = "".join(symbolic_parts) + " =" + "".join(value_parts) + f" = {total}"

    print(equation_str)
    
    return total

# Main execution block
if __name__ == '__main__':
    # You can change these values to solve for different parameters.
    # The problem specifies N, K, M are positive integers and M*(K-1) < N.
    N = 20
    K = 5
    M = 4
    
    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}")
    print("The formula used is: Sum_{p=0 to K-1} [(-1)^p * C(K-1, p) * C(N - p*M, K)]")
    print("Equation:")
    
    # Calculate and print the result
    final_result = calculate_sequences(N, K, M)
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"\nTotal number of possible sequences: {final_result}")
    # The submission format requires the answer at the end.
    # So I will not print it here but put it at the end of the response.