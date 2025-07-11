import math

def calculate_mz1(n):
    """
    Calculates M_z(1, n) based on the derived formula.
    
    The formula is:
    M_z(1, n) = (-2/(pi*n))^n * (4n)!/n! * sum
    where sum = sum_{k=0}^{n+1} C(n+1, k) * (-1)^k / (3n-1+k)!
    """
    if n <= 0:
        return float('inf')

    # Prefactor calculation using large integer arithmetic
    try:
        prefactor_log = n * (math.log(2) - math.log(math.pi) - math.log(n)) \
                      + math.lgamma(4 * n + 1) - math.lgamma(n + 1)
        
        # Determine the sign
        sign = -1 if n % 2 != 0 else 1
        
        prefactor = sign * math.exp(prefactor_log)
    except ValueError:
        return float('inf') # Handles cases like log(0) or lgamma of non-positive

    # Sum calculation
    sum_val = 0
    for k in range(n + 2):
        # Calculate term k of the sum
        # C(n+1, k) * (-1)^k / (3n-1+k)!
        try:
            log_term = math.lgamma(n + 2) - math.lgamma(k + 1) - math.lgamma(n + 2 - k) \
                       - math.lgamma(3 * n + k)
            
            term_sign = -1 if k % 2 != 0 else 1
            
            term = term_sign * math.exp(log_term)
            sum_val += term
        except ValueError:
            # This can happen if 3n+k is not positive, but it won't for n>=1
            continue

    return prefactor * sum_val

def main():
    """
    Finds the n that minimizes M_z(1, n) and prints the calculation.
    """
    min_mz = float('inf')
    n_min = -1
    mz_values = []

    # Calculate M_z(1, n) for n from 1 to 20
    for n in range(1, 21):
        mz = calculate_mz1(n)
        mz_values.append((n, mz))
        if mz < min_mz:
            min_mz = mz
            n_min = n

    print(f"Searching for the minimum magnetization M_z(1) for n = 1, 2, ..., 20.")
    for n, mz in mz_values:
        print(f"For n = {n}, M_z(1) = {mz:.6f}")

    print(f"\nThe minimum magnetization occurs at n = {n_min}.")
    print(f"The minimum value is M_z(1, {n_min}) = {min_mz:.6f}\n")

    # Output the detailed calculation for n_min
    n = n_min
    pi = math.pi
    
    print(f"The detailed calculation for n = {n}:")
    print(f"M_z(1, {n}) = (-2 / (pi * {n}))^{n} * ({4*n}! / {n}!) * Sum")
    
    # Calculate prefactor parts
    p1 = f"(-2 / ({pi:.4f} * {n}))^{n}"
    p2 = f"({4*n}! / {n}!)"
    
    sum_str = []
    sum_val_num = 0
    sum_val_den = math.factorial(3 * n - 1 + n + 1) # Common denominator for display is complex
    
    sum_terms_str = []
    
    sum_val = 0
    for k in range(n + 2):
        sign = (-1)**k
        comb_val = math.comb(n + 1, k)
        fact_val = math.factorial(3 * n - 1 + k)
        term = sign * comb_val / fact_val
        sum_val += term
        
        sign_str = "-" if sign < 0 else "+"
        if k==0: sign_str = ""
            
        sum_terms_str.append(f" {sign_str} C({n+1},{k}) / ({3*n-1+k})!")

    print("where Sum = " + "".join(sum_terms_str))

    # Calculate the final result using high precision math
    prefactor = ((-2 / (pi * n))**n) * (math.factorial(4*n) / math.factorial(n))
    final_result = prefactor * sum_val

    print(f"\nEvaluating the expression:")
    print(f"Prefactor = {p1} * {p2} = {prefactor:.6f}")
    print(f"Sum = {sum_val:.6e}")
    print(f"M_z(1, {n}) = {prefactor:.6f} * {sum_val:.6e} = {final_result:.6f}")

if __name__ == "__main__":
    main()
