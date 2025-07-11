import math

def solve_simplices_count():
    """
    Calculates and prints the number of n-simplices for N=200, k=13, and n<=5.
    """
    # Given parameters from the problem
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N} and k={k}.\n")
    
    # Store results for the final answer
    results = []

    # Loop for n from 0 to 5
    for n in range(6):
        # The number of n-simplices is given by C(N-k+n+1, n+1)
        a = N - k + n + 1
        b = n + 1

        # Calculate the result using math.comb for accuracy and efficiency
        result = math.comb(a, b)
        results.append(result)

        print(f"For n = {n}:")
        
        # Build the numerator and denominator strings for the expanded formula
        # The numerator is a * (a-1) * ... * (a-b+1)
        numerator_vals = [str(i) for i in range(a, a - b, -1)]
        numerator_str = " * ".join(numerator_vals)
        
        # The denominator is b! = b * (b-1) * ... * 1
        denominator_vals = [str(i) for i in range(b, 0, -1)]
        denominator_str = " * ".join(denominator_vals)

        print(f"The number is C({a}, {b}) = ({numerator_str}) / ({denominator_str}) = {result}\n")
    
    return results

if __name__ == '__main__':
    solve_simplices_count()
