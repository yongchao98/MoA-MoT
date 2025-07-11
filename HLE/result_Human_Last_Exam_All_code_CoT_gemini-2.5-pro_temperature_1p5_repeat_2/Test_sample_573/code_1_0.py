import math

def solve():
    """
    Calculates the number of n-simplices for n<=5, N=200, and k=13.
    """
    N = 200
    k = 13
    
    # The number of objects in the relevant subcategory is p
    p = N - k + 1

    print(f"Given N = {N} and k = {k}.")
    print(f"The number of n-simplices is given by the formula C(p + n, n + 1), where p = N - k + 1 = {p}.\n")

    results = []
    for n in range(6):
        # The number of items to choose is m = n + 1
        m = n + 1
        
        # The formula for combinations with repetition is C(p + m - 1, m)
        # C(p + (n+1) - 1, n+1) = C(p+n, n+1)
        
        # In our case, this is C(188 + n, n+1)
        upper_val = p + n
        lower_val = m
        
        # Calculate the number of simplices using math.comb
        num_simplices = math.comb(upper_val, lower_val)
        
        print(f"For n = {n}:")
        equation = f"  The number of {n}-simplices is C({N} - {k} + {n} + 1, {n} + 1) = C({upper_val}, {lower_val}) = {num_simplices}"
        print(equation)
        results.append(num_simplices)

    # Returning the final numerical answer as requested by the wrapper.
    # The prompt implies a final single answer is expected in the <<<>>> format.
    # The problem asks "How many n-simplices ... for n<=5", suggesting a list of answers.
    # The values will be formatted into the answer string.
    
    answer_string = ", ".join([f"n={i}: {res}" for i, res in enumerate(results)])
    # print(f"\n<<< {answer_string} >>>")

solve()
