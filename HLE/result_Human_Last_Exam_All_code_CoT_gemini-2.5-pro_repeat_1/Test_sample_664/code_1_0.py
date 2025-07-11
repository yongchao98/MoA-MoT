import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def double_factorial(n):
    """Calculates the double factorial n!!"""
    if n < -1:
        raise ValueError("Double factorial is not defined for numbers less than -1")
    if n == -1 or n == 0:
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

def count_involutions(n):
    """
    Counts the number of involutions on n elements.
    An involution is a permutation that is its own inverse.
    The number of involutions I(n) is the sum of C(n, 2k) * (2k-1)!! for k=0 to n//2.
    This corresponds to choosing 2k elements to be in k pairs (2-cycles),
    and the remaining n-2k elements are fixed points (1-cycles).
    """
    total_involutions = 0
    terms = []
    
    # k is the number of 2-cycles (transpositions)
    for k in range(n // 2 + 1):
        num_elements_in_2_cycles = 2 * k
        
        # Choose 2k elements out of n
        ways_to_choose = combinations(n, num_elements_in_2_cycles)
        
        # Pair them up. The number of ways to form k pairs from 2k elements is (2k-1)!!
        ways_to_pair = double_factorial(num_elements_in_2_cycles - 1)
        
        term = ways_to_choose * ways_to_pair
        terms.append(term)
        total_involutions += term

    # Format the output equation
    equation_parts = []
    for k in range(n // 2 + 1):
        num_elements = 2 * k
        c_nk = combinations(n, num_elements)
        df = double_factorial(num_elements - 1)
        # We can simplify the output to just show the value of each term
        # equation_parts.append(f"C({n},{num_elements}) * {num_elements-1}!! = {c_nk} * {df} = {terms[k]}")
        equation_parts.append(str(terms[k]))
    
    print(f"The number of configurations is the number of involutions on 8 elements, I(8).")
    print(f"This is calculated by summing the ways to have k pairs of swapped chips (2-cycles) for k=0 to 4.")
    print("I(8) = (ways for 0 pairs) + (ways for 1 pair) + (ways for 2 pairs) + (ways for 3 pairs) + (ways for 4 pairs)")
    
    final_equation = " + ".join(equation_parts)
    print(f"I(8) = {final_equation}")
    print(f"The total number of possible configurations is {total_involutions}.")


if __name__ == "__main__":
    count_involutions(8)
<<<764>>>