import math

def calculate_bound(n, s):
    """
    Calculates the bound m <= sum_{i=0 to s} C(n-1, i).

    Args:
        n (int): The size of the base set [n].
        s (int): The number of integers in the set L.
    
    Returns:
        A tuple containing the bound and the list of terms in the sum.
    """
    if n <= 0 or s < 0 or s >= n:
        raise ValueError("Inputs must satisfy n > 0, 0 <= s < n.")
    
    # The bound is on subsets of [n-1]
    N = n - 1
    
    total_sum = 0
    terms = []
    for i in range(s + 1):
        # Calculate C(N, i)
        term = math.comb(N, i)
        terms.append(term)
        total_sum += term
        
    return total_sum, terms

def main():
    """
    Main function to demonstrate the bound calculation.
    """
    # Example values from the problem context
    n = 10
    s = 3
    
    try:
        bound, terms = calculate_bound(n, s)
        
        # Format the output equation
        terms_str = " + ".join(map(str, terms))
        
        print(f"For n = {n} and s = {s}, the bound from the theorem in part (b) is:")
        print(f"m <= sum_{{i=0}}^{s} C(n-1, i)")
        print(f"m <= C({n-1}, 0) + ... + C({n-1}, {s})")
        print("Plugging in the numbers:")
        print(f"m <= {terms_str}")
        print(f"m <= {bound}")
        
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
