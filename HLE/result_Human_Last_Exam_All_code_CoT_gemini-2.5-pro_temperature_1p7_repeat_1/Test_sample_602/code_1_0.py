import math

def calculate_l(n):
    """
    Calculates the value of l(n) for a given integer n >= 5.

    The function first prints the steps of the calculation using the derived
    analytical formula for l(n), and then prints the final numerical result.

    The formula is: l(n) = 2 + 2/n^2 - (2*(2n-1)*sqrt(n^2-n+1))/n^2
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    print(f"Calculating l(n) for n = {n}")
    
    # Define intermediate variables for clarity in printing
    n_squared = n * n
    term_1 = 2
    term_2_num = 2
    term_2 = term_2_num / n_squared
    
    sqrt_term_val = math.sqrt(n_squared - n + 1)
    term_3_num_factor1 = 2 * (2 * n - 1)
    term_3_num = term_3_num_factor1 * sqrt_term_val
    term_3 = term_3_num / n_squared
    
    result = term_1 + term_2 - term_3

    # Printing the breakdown of the final equation
    print("The final equation for l(n) is: l(n) = 2 + 2/n^2 - (2 * (2n - 1) * sqrt(n^2 - n + 1)) / n^2")
    print(f"Substituting n = {n}:")
    print(f"l({n}) = 2 + 2/({n}^2) - (2 * (2*{n} - 1) * sqrt({n}^2 - {n} + 1)) / ({n}^2)")
    print(f"l({n}) = {term_1} + {term_2_num}/{n_squared} - (2 * {2*n-1} * sqrt({n_squared - n + 1})) / {n_squared}")
    print(f"l({n}) = {term_1} + {term_2} - ({term_3_num_factor1} * {sqrt_term_val}) / {n_squared}")
    print(f"l({n}) = {term_1 + term_2} - {term_3_num} / {n_squared}")
    print(f"l({n}) = {term_1 + term_2} - {term_3}")
    print(f"Final Result: l({n}) = {result}")

if __name__ == '__main__':
    # Example usage for n=5
    calculate_l(5)
    print("\n" + "="*20 + "\n")
    # Example usage for n=10
    calculate_l(10)