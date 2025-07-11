import math

def double_factorial(n):
    """Computes the double factorial n!!."""
    if n < 0:
        return 0
    if n == 0 or n == 1:
        return 1
    res = 1
    for i in range(n, 0, -2):
        res *= i
    return res

def compute_chi_mod_k(k):
    """
    Computes the reduced Euler characteristic of Delta_k modulo k.
    
    The formula is: (-1)^((k+1)/2) * (k-2)!! mod k
    """
    if k < 3 or k % 2 == 0:
        # The problem states k is a prime >= 3, which means k is odd.
        # This is a check for valid input based on the problem statement.
        print("k must be an odd prime >= 3.")
        return

    # Step 1: Calculate the sign
    power_of_neg_one = (k + 1) // 2
    sign = -1 if power_of_neg_one % 2 != 0 else 1

    # Step 2: Calculate (k-2)!! mod k
    # We can do this iteratively to show the steps.
    double_fact_val = 1
    double_fact_terms = []
    if k > 2:
        for i in range(k - 2, 0, -2):
            double_fact_val = (double_fact_val * i)
            double_fact_terms.append(str(i))
    
    double_fact_mod_k = double_fact_val % k

    # Step 3: Combine and find the final result
    result = (sign * double_fact_mod_k) % k

    # Print the equation
    df_str = " * ".join(double_fact_terms)
    if not df_str:
        df_str = "1"
        
    print(f"For k = {k}:")
    print(f"  Result = (-1)^(({k}+1)/2) * ({k}-2)!! mod {k}")
    print(f"         = ({sign}) * ({df_str}) mod {k}")
    print(f"         = ({sign}) * ({double_fact_val}) mod {k}")
    print(f"         = ({sign}) * ({double_fact_mod_k}) mod {k}")
    print(f"         = {sign * double_fact_mod_k} mod {k}")
    print(f"         = {result}")
    
    # Return the final integer result for the "answer" block
    return result

# You can change the value of k here to test other primes
k = 7
final_answer = compute_chi_mod_k(k)
# The final answer format, as requested
print(f"<<<{final_answer}>>>")
