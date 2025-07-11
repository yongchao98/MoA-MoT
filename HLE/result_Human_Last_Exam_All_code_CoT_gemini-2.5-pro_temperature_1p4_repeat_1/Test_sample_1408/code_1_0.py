import math

def solve_overhang_puzzle():
    """
    This function calculates the integers a, b, c for the maximal 3-cube overhang problem.
    The maximal overhang is found to be 5*sqrt(2)/6 using a 'diamond' configuration.
    We need to find integers a,b,c such that (a + sqrt(b)) / (1 + c) = 5*sqrt(2)/6,
    with c being the minimal non-negative integer.
    """

    # The target value is 5*sqrt(2)/6
    # From derivation, for (a+sqrt(b))/(1+c) to equal the target, we must have a=0.
    # This leads to b = (25 * (1+c)^2) / 18.
    # We need to find the smallest non-negative integer c for which b is an integer.

    c = 0
    while True:
        # Check if (1+c)^2 is divisible by 18, so that b is an integer.
        # This simplifies to checking if (1+c) is a multiple of 6.
        if (1 + c) % 6 == 0:
            a = 0
            b_numerator = 25 * (1 + c)**2
            b_denominator = 18
            if b_numerator % b_denominator == 0:
                b = b_numerator // b_denominator

                # Check if sqrt(b) is non-integer
                sqrt_b = math.isqrt(b)
                if sqrt_b * sqrt_b != b:
                    # We have found the minimal c and corresponding a, b.
                    print("The maximal overhang for three cubes can be expressed in the form (a + sqrt(b))/(1 + c).")
                    print("Based on the optimal 'diamond' configuration, the overhang is 5*sqrt(2)/6.")
                    print("Solving for the parameters with the constraint that c is minimal, we get:")
                    print(f"a = {a}")
                    print(f"b = {b}")
                    print(f"c = {c}")
                    print("\nThe equation is:")
                    print(f"({a} + sqrt({b})) / (1 + {c}) = {a} + sqrt({b}) / {1+c}")
                    return a, b, c
        c += 1

# Execute the solver
a, b, c = solve_overhang_puzzle()
# The final output format required
# <<<a b c>>>
# print(f'<<<{a} {b} {c}>>>')
