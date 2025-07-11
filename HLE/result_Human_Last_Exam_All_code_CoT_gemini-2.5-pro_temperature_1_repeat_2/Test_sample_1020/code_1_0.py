import math

def combinations(n, k):
    """
    Computes the binomial coefficient C(n, k), also known as "n choose k".
    """
    if k < 0 or k > n:
        return 0
    # math.comb is efficient and handles large numbers.
    return math.comb(n, k)

def solve_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: The number of available indices.
        d: The degree of the polynomial (number of variables in each monomial).
    
    Returns:
        A string representing the equation for the complexity.
    """
    if not (2 <= d <= n):
        raise ValueError("The condition 2 <= d <= n must be met.")

    terms = [2]
    total_complexity = 2
    
    for i in range(1, d):
        term = combinations(n, i)
        terms.append(term)
        total_complexity += term
        
    equation_parts = [str(t) for t in terms]
    equation = " + ".join(equation_parts) + f" = {total_complexity}"
    
    return equation

def main():
    """
    Main function to demonstrate the solution with example values.
    """
    # Example values for n and d, where d is odd.
    n = 6
    d = 5
    
    print(f"For n = {n} and d = {d}:")
    
    # Calculate the smallest complexity
    complexity_equation = solve_complexity(n, d)
    
    # Print the result in the specified format
    print(f"The smallest complexity is given by the sum of the dimensions of the intermediate matrices plus 2.")
    print(f"The complexity is calculated as: 2 + C({n},1) + C({n},2) + ... + C({n},{d-1})")
    print("The final equation is:")
    print(complexity_equation)

if __name__ == "__main__":
    main()
