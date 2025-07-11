import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k), also known as "n choose k".
    """
    if k < 0 or k > n:
        return 0
    # This is an efficient way to calculate n! / (k! * (n-k)!)
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_a_n(n):
    """
    Calculates a(n), the maximal number of prime implicants of a Boolean function of n variables.
    It prints the step-by-step calculation.
    """
    total_sum = 0
    terms = []
    
    print(f"To find a({n}), we use the formula: a(n) = sum for k=0 to {n} of [C(n, k) * C(k, floor(k/2))]")
    print("-" * 20)

    for k in range(n + 1):
        # Calculate the two combination terms
        c_nk = combinations(n, k)
        k_floor_div_2 = k // 2
        c_k_floor = combinations(k, k_floor_div_2)
        
        # Calculate the product for the current k
        term_value = c_nk * c_k_floor
        
        # Print the breakdown for the current term
        print(f"Term for k={k}: C({n}, {k}) * C({k}, {k_floor_div_2}) = {c_nk} * {c_k_floor} = {term_value}")
        
        # Add to the total sum and store the term for the final equation
        total_sum += term_value
        terms.append(str(term_value))

    # Construct and print the final equation
    final_equation = " + ".join(terms)
    print("\n" + "-" * 20)
    print(f"So, a({n}) = {final_equation}")
    print(f"a({n}) = {total_sum}")


# The user wants to know a(4)
n_value = 4
calculate_a_n(n_value)
