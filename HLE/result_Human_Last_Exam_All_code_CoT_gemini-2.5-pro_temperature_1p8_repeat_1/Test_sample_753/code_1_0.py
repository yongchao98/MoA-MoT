import math
from decimal import Decimal, getcontext

def solve():
    """
    Solves the problem by calculating the first 40 digits of the cardinality of Theta^{-1}(lambda) for m=3.
    """
    
    # Set the value of m as given in the problem
    m = 3
    
    # Step 1: Determine n and the partition lambda for m=3
    print(f"For m = {m}, we first calculate n using the formula n = sum_{k=1 to m} k*(m+1-k):")
    n = 0
    n_calc_str_parts = []
    partition_desc = {}
    for k in range(1, m + 1):
        term = k * (m + 1 - k)
        n += term
        n_calc_str_parts.append(f"{k}*({m}+1-{k})")
        
        # Determine the parts of the partition lambda
        part_value = m + 1 - k
        part_multiplicity = k
        if part_value in partition_desc:
            partition_desc[part_value] += part_multiplicity
        else:
            partition_desc[part_value] = part_multiplicity
    
    print(f"n = {' + '.join(n_calc_str_parts)} = {sum(k*(m+1-k) for k in range(1,m+1))}")
    print(f"So, we are working with the symmetric group Sigma_{n}, where n = {n}.")
    
    lambda_str = ", ".join([f"{key}^{val}" for key, val in sorted(partition_desc.items(), reverse=True)])
    print(f"\nThe partition lambda is given by (m^1, (m-1)^2, ..., 1^m), which for m={m} is ({lambda_str}).")

    # Step 2: Calculate the size of the conjugacy class C_lambda
    n_factorial = Decimal(math.factorial(n))
    
    # For m=3, the partition is (3^1, 2^2, 1^3).
    # This corresponds to cycle counts c_3=1, c_2=2, c_1=3.
    # From the general formula, we get the same counts:
    cycle_counts = { m+1-k : k for k in range(1,m+1)}
    
    # Formula for size of conjugacy class: n! / (product of k^c_k * c_k!)
    denom = Decimal(1)
    denom_str_parts = []
    for k, c_k in sorted(cycle_counts.items(), reverse=True):
        term = (Decimal(k)**c_k) * Decimal(math.factorial(c_k))
        denom *= term
        denom_str_parts.append(f"{k}^{c_k}*{c_k}!")
    
    c_lambda_size = n_factorial / denom
    
    print(f"\nThe size of the conjugacy class C_lambda is calculated as:")
    print(f"|C_lambda| = {n}! / ( {' * '.join(denom_str_parts)} )")
    print(f"|C_lambda| = {n_factorial} / {denom} = {c_lambda_size.to_integral_value()}")
    
    # Step 3: High-precision calculation for the cardinality
    # Based on our analysis, the cardinality of Theta^{-1}(lambda) is |C_lambda| * (n!)^(n!).
    cardinality_coeff = c_lambda_size
    cardinality_base = n_factorial
    cardinality_exp = n_factorial

    print("\nThe cardinality of Theta^-1(lambda) is given by the final equation:")
    final_eq_str = f"|C_lambda| * (n!)^(n!) = {cardinality_coeff.to_integral_value()} * {cardinality_base.to_integral_value()}^{cardinality_base.to_integral_value()}"
    print(final_eq_str)

    # To find the first 40 digits, we compute the base-10 logarithm of this number.
    # log10(Result) = log10(coeff) + exp * log10(base).
    # The fractional part 'f' of this logarithm gives the leading digits via 10^f.
    getcontext().prec = 100 # Set precision high enough for the calculation.

    log10_coeff = cardinality_coeff.log10()
    log10_base = cardinality_base.log10()
    
    log10_result = log10_coeff + cardinality_exp * log10_base
    
    fractional_part = log10_result - log10_result.to_integral_value(rounding='ROUND_FLOOR')
    
    # The leading digits are given by 10^fractional_part
    first_digits_val = Decimal(10).pow(fractional_part)

    # Format the output to get the first 40 digits.
    # We generate a string with enough precision and then truncate.
    first_digits_str = ('%.45f' % first_digits_val).replace('.', '')[:40]

    print("\nThe first 40 digits of this cardinality are:")
    print(first_digits_str)

# Execute the solver function
solve()