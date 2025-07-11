import math

def solve():
    """
    Solves the problem by calculating the required limit value.
    
    The plan is as follows:
    1.  The problem's conditions imply that for any solution pair (a,b), the expression
        k = (a^2 + b^2 + 5a + 5b + 1) / (ab) must be an integer.
    2.  Using a proof technique called Vieta jumping, we can demonstrate that this integer k
        is constant for all solutions, and its value is 13.
    3.  This means all solutions (a,b) satisfy the Diophantine equation:
        a^2 + b^2 + 5a + 5b + 1 = 13ab.
    4.  The solutions to this equation form an infinite sequence. The growth rate of this
        sequence is determined by the larger root of the characteristic equation r^2 - 13r + 1 = 0.
    5.  The number of solutions F(N) for a, b <= N grows proportionally to ln(N). The limit
        lim_{N->inf} F(N)/ln(N) can be calculated based on this growth rate.
    6.  The limit, L, is found to be 2 / ln(alpha), where alpha is the larger root of the
        characteristic equation.
    7.  The final result is the integer part of 10^4 * L.
    """
    
    # Numbers from the characteristic equation r^2 - 13r + 1 = 0
    # which governs the growth of the solution sequence.
    char_eq_coeff_b = 13
    char_eq_coeff_c = 1
    
    # Calculate the discriminant: b^2 - 4ac
    discriminant = char_eq_coeff_b**2 - 4 * 1 * char_eq_coeff_c
    
    # The growth rate 'alpha' is the larger root of the characteristic equation
    alpha = (char_eq_coeff_b + math.sqrt(discriminant)) / 2
    
    # The limit L = lim_{N->inf} F(N)/ln(N) is given by 2 / ln(alpha)
    limit_numerator = 2
    limit_L = limit_numerator / math.log(alpha)
    
    # The problem asks for the integer part of 10^4 * L
    scaling_factor = 10000
    final_value = scaling_factor * limit_L
    
    integer_part = int(final_value)

    # Output the steps of the calculation with the numbers involved as requested.
    print(f"The growth rate of the solutions is determined by the root of the equation r^2 - {char_eq_coeff_b}r + {char_eq_coeff_c} = 0.")
    print(f"The discriminant is {char_eq_coeff_b}^2 - 4*1*{char_eq_coeff_c} = {discriminant}.")
    print(f"The growth rate alpha = ({char_eq_coeff_b} + sqrt({discriminant})) / 2.")
    print(f"The limit L is {limit_numerator} / ln(alpha).")
    print(f"The final expression to calculate is: floor({scaling_factor} * {limit_numerator} / ln(({char_eq_coeff_b} + sqrt({discriminant})) / 2))")
    print(f"The computed value is: {final_value}")
    print(f"The integer part is: {integer_part}")

solve()