import math

def solve():
    """
    Calculates the final value based on the derivation.
    """
    # p is the Mersenne prime 2^127 - 1
    p = (2**127) - 1

    # The problem reduces to computing -C^2 mod p.
    # C is the product of the modular multinomial coefficients for the base-p digits.
    # The repeating factors are for digit combinations (1,4,1), (3,2,3), and (4,2,4).

    # Factor for digits (1, 4, 1): (1+4+1)! / (1! * 4! * 1!)
    f1 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))

    # Factor for digits (3, 2, 3): (3+2+3)! / (3! * 2! * 3!)
    f2 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))

    # Factor for digits (4, 2, 4): (4+2+4)! / (4! * 2! * 4!)
    f3 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

    # C is the product of these factors
    C = f1 * f2 * f3

    # The result is p - C^2
    C_squared = C * C
    result = p - C_squared

    print("The final calculation is based on the equation: result = p - C^2")
    print(f"The value of p is: {p}")
    print(f"The value of C = {f1} * {f2} * {f3} = {C}")
    print(f"The value of C^2 is: {C_squared}")
    print(f"The final result is {p} - {C_squared} which equals:")
    print(result)

solve()