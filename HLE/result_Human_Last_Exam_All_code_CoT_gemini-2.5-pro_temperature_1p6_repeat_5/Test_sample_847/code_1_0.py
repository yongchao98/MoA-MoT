import math

def main():
    """
    This script calculates the minimal cost per bit of information to determine
    the permutation of an array under the given questioning rules.
    """
    
    # The optimal strategy for Type 2 questions requires balancing the cost per
    # bit of information for "yes" and "no" answers. This leads to the equation:
    # 1*p^3 + 1*p - 1 = 0
    # where p is the optimal probability of a "yes" answer.
    # We output the coefficients of this equation as requested.
    
    coeff_p3 = 1
    coeff_p1 = 1
    const_term = -1
    
    print("The optimal probability 'p' is the root of the equation:")
    print(f"({coeff_p3}) * p^3 + ({coeff_p1}) * p + ({const_term}) = 0")
    print("-" * 30)

    # We find the root p0 of p^3 + p - 1 = 0 using the bisection method.
    def p_equation(p):
        return p**3 + p - 1

    a, b = 0.0, 1.0
    for _ in range(100):  # 100 iterations provide high precision
        mid = (a + b) / 2
        if p_equation(mid) < 0:
            a = mid
        else:
            b = mid
    p0 = (a + b) / 2

    # Calculate the cost per bit for the two types of questions.
    # Cost per bit for Type 1 is fixed at 2.
    cost_per_bit_type1 = 2.0
    
    # Cost per bit for Type 2 is -1 / log2(p0).
    cost_per_bit_type2 = -1 / math.log2(p0)

    # The minimal cost per bit is the minimum of the two.
    min_cost_per_bit = min(cost_per_bit_type1, cost_per_bit_type2)
    
    print(f"Cost per bit using Type 1 questions: {cost_per_bit_type1:.3f}")
    print(f"Optimal cost per bit using Type 2 questions: {cost_per_bit_type2:.3f}")
    print("-" * 30)
    print("The minimal cost per bit of information for large n is the minimum of these two values.")
    print(f"Minimal cost per bit = {min_cost_per_bit:.3f}")

if __name__ == '__main__':
    main()
