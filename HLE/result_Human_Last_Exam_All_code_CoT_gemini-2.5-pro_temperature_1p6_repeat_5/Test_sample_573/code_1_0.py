import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for N=200, k=13, and n from 0 to 5.
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N} and k={k}.")
    print("The formula for the number of n-simplices is C(N - k + n + 1, n + 1).\n")
    
    results = {}

    for n in range(6):
        upper = N - k + n + 1
        lower = n + 1
        
        result = math.comb(upper, lower)
        results[n] = result
        
        print(f"For n = {n}:")
        print(f"Number of {n}-simplices = C({N} - {k} + {n} + 1, {n} + 1)")
        
        if lower > 1:
            numerator_parts = [str(i) for i in range(upper, upper - lower, -1)]
            denominator_parts = [str(i) for i in range(lower, 0, -1)]
            
            numerator_str = " * ".join(numerator_parts)
            denominator_str = " * ".join(denominator_parts)
            
            print(f"                    = C({upper}, {lower}) = ({numerator_str}) / ({denominator_str})")
        else:
            print(f"                    = C({upper}, {lower})")

        print(f"                    = {result}\n")
        
    return results

if __name__ == "__main__":
    final_results = calculate_simplices()
    # The final answer format is specified in the problem.
    # To conform to this, we will prepare the dictionary string for the final output.
    # print(f"\nFinal answer dictionary: {final_results}")