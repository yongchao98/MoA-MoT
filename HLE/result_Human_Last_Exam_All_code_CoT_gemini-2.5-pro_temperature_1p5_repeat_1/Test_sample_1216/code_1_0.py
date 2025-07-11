import math

def solve_augustus_records():
    """
    Solves the Augustus records problem by finding integer solutions
    to the system of equations derived from the problem's constraints.
    """
    # The sum D + S + F must equal 400, as derived from the problem statement:
    # Total readable records = 720
    # ICDF records = 720 / 3 = 240
    # Caesar-only records = 80
    # Remainder for D, S, F = 720 - 240 - 80 = 400
    total_for_dsf = 400
    
    # We search for d = sqrt(D) and s = sqrt(S).
    # d and s must be integers. s^2 < 400 means s < 20.
    s_sol, d_sol = -1, -1

    # We are looking for an integer solution for 'd' in the quadratic equation:
    # d^2 + d + (s^2 + 12s - 400) = 0
    # For 'd' to be an integer, the discriminant must be a perfect square.
    # Discriminant = 1^2 - 4*1*(s^2 + 12s - 400) = 1601 - 4s^2 - 48s
    
    for s_candidate in range(1, int(math.sqrt(total_for_dsf))):
        discriminant = 1601 - 4 * s_candidate**2 - 48 * s_candidate
        if discriminant >= 0:
            sqrt_discriminant = math.isqrt(discriminant)
            if sqrt_discriminant**2 == discriminant:
                # Check if this leads to a positive integer solution for d
                # d = (-b + sqrt(discriminant)) / 2a
                numerator = -1 + sqrt_discriminant
                if numerator > 0 and numerator % 2 == 0:
                    d_candidate = numerator // 2
                    # Verify our findings
                    if d_candidate**2 + s_candidate**2 + d_candidate + 12 * s_candidate == total_for_dsf:
                        s_sol = s_candidate
                        d_sol = d_candidate
                        break

    if s_sol == -1:
        print("Could not find an integer solution.")
        return

    # Calculate the final values based on the solution
    lost_records_l = 4 * s_sol
    single_variant_s = s_sol**2
    dual_named_d = d_sol**2
    full_title_f = d_sol + 12 * s_sol

    # The denominator is the number of distinct patterns containing the name "Augustus"
    # 1. "Augustus" and "Caesar" (D)
    # 2. "Augustus" alone (from S)
    # 3. "Imperator Caesar Augustus" (F)
    denominator = 3
    
    # Calculate the product for the numerator
    product = lost_records_l * dual_named_d * single_variant_s * full_title_f
    
    # Final calculation steps
    result = product / denominator
    final_answer = math.floor(result / 1000)

    # Print the required outputs
    print(f"Based on the problem's constraints, the calculated values are:")
    print(f"a) Lost records (L): {lost_records_l}")
    print(f"b) Dual-named documents (D): {dual_named_d}")
    print(f"c) Single-variant documents (S): {single_variant_s}")
    print(f"d) Full imperial title documents (F): {full_title_f}\n")
    
    print("Final equation steps:")
    print(f"Product = {lost_records_l} * {dual_named_d} * {single_variant_s} * {full_title_f} = {product}")
    print(f"Denominator (Number of distinct Augustus naming patterns) = {denominator}")
    print(f"Result = {product} / {denominator} = {int(result)}")
    print(f"Final Answer = floor({int(result)} / 1000) = {final_answer}")

    print(f"\n<<<{final_answer}>>>")

solve_augustus_records()