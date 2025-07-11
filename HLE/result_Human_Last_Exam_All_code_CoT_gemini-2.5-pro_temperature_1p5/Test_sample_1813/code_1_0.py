import math
from decimal import Decimal, getcontext

def compute_cfrac_for_markov_number(r, s):
    """
    Computes and prints the continued fraction for the generalized Markov number m_{r/s}.
    """
    # 1. Calculate the value of the Markov number
    # m_{r/s} = (sqrt(r^2 + 4s^2))/s
    d_inside_sqrt = r**2 + 4 * s**2
    
    # We look for a simplified form a*sqrt(d)/c
    # We can pull out perfect square factors from d_inside_sqrt
    # Let d_inside_sqrt = factor^2 * d
    factor = int(math.sqrt(d_inside_sqrt))
    d = d_inside_sqrt
    simplified_b = 1
    while factor > 1:
        if d_inside_sqrt % (factor**2) == 0:
            d = d_inside_sqrt // (factor**2)
            simplified_b = factor
            break
        factor -= 1

    # Our number is simplified_b * sqrt(d) / s
    # We will compute the continued fraction of (a_start + b_start * sqrt(d)) / c_start
    a_start = 0
    b_start = simplified_b
    c_start = s

    # 2. Print the derivation
    print("Step 1: Calculate the value of the generalized Markov number m_{r/s}")
    print("The formula is: m_{r/s} = (sqrt(r^2 + 4*s^2)) / s")
    print(f"For r={r} and s={s}, we get:")
    print(f"m_{r}/{s} = (sqrt({r}^2 + 4*{s}^2)) / {s}")
    print(f"      = (sqrt({r**2} + 4*{s**2})) / {s}")
    print(f"      = (sqrt({r**2 + 4*s**2})) / {s}")
    if simplified_b > 1:
        print(f"      = (sqrt({simplified_b**2} * {d})) / {s}")
        print(f"      = ({simplified_b}*sqrt({d})) / {s}")
    print("\nStep 2: Compute the continued fraction of this number.\n")

    # 3. Compute the continued fraction using integer arithmetic
    getcontext().prec = 100 # Set precision for floating point estimation

    a_curr, b_curr, c_curr = a_start, b_start, c_start
    
    coeffs = []
    history = {} # To detect periodicity
    
    for _ in range(100): # Safety limit for iterations
        state = (a_curr, b_curr, c_curr)
        if state in history:
            period_start_index = history[state]
            non_periodic_part = coeffs[:period_start_index]
            periodic_part = coeffs[period_start_index:]
            break
        
        history[state] = len(coeffs)

        # Estimate and find the integer part
        val_estimate = (Decimal(a_curr) + Decimal(b_curr) * Decimal(d).sqrt()) / Decimal(c_curr)
        integer_part = int(val_estimate)
        coeffs.append(integer_part)

        # Compute the next remainder: 1 / (x - k)
        # x - k = (a + b*sqrt(d))/c - k = (a - k*c + b*sqrt(d))/c
        # 1/(x-k) = c / (a-k*c + b*sqrt(d))
        # Rationalize: c * (a-k*c - b*sqrt(d)) / ((a-k*c)^2 - b^2*d)
        
        a_temp = a_curr - integer_part * c_curr
        
        # New components for the next iteration
        a_next = c_curr * a_temp
        b_next = -c_curr * b_curr
        c_next = a_temp**2 - b_curr**2 * d
        
        # Simplify by dividing by the greatest common divisor
        common_divisor = math.gcd(math.gcd(a_next, b_next), c_next)
        
        a_curr = a_next // common_divisor
        b_curr = b_next // common_divisor
        c_curr = c_next // common_divisor
        
        # Conventionally, the denominator should be positive
        if c_curr < 0:
            a_curr = -a_curr
            b_curr = -b_curr
            c_curr = -c_curr
    else:
        # This part should not be reached for quadratic irrationals
        non_periodic_part = coeffs
        periodic_part = []

    # 4. Print the final result
    print("The continued fraction expansion for m_{4/7} is:")
    non_periodic_str = ", ".join(map(str, non_periodic_part))
    periodic_str = ", ".join(map(str, periodic_part))

    final_equation = f"m_{r}/{s} = [{non_periodic_str}"
    if periodic_part:
        final_equation += f"; ({periodic_str})]"
    else:
        final_equation += "]"
        
    print(final_equation)
    # The user asked to output each number, let's list them explicitly too.
    print("\nThe coefficients are:")
    print(f"a_0 = {non_periodic_part[0]}")
    for i, p in enumerate(periodic_part, 1):
        print(f"a_{i} = {p}")

if __name__ == '__main__':
    # We want to compute the continued fraction for m_{4/7}
    compute_cfrac_for_markov_number(r=4, s=7)