import math

def solve_a4():
    """
    Calculates a(4), the maximal number of prime implicants of a Boolean function of 4 variables.

    The formula for a(n) is:
    a(n) = sum_{k=0 to n} C(n, k) * C(k, floor(k/2))
    where C(n, k) is the binomial coefficient "n choose k".

    We will compute this for n = 4.
    """
    n = 4
    terms = []
    total_sum = 0

    # Loop from k=0 to n to calculate each term in the summation
    for k in range(n + 1):
        # Calculate C(n, k)
        c_n_k = math.comb(n, k)
        
        # Calculate C(k, floor(k/2))
        floor_k_div_2 = k // 2
        c_k_floor = math.comb(k, floor_k_div_2)
        
        # The term is the product of the two binomial coefficients
        term = c_n_k * c_k_floor
        terms.append(term)
        total_sum += term

    # Construct the final equation string as requested
    equation_str = " + ".join(map(str, terms))
    final_output = f"{equation_str} = {total_sum}"

    print("To find a(4), we calculate the sum of C(4,k) * C(k, floor(k/2)) for k from 0 to 4.")
    print("The final equation is the sum of these terms:")
    print(final_output)

solve_a4()