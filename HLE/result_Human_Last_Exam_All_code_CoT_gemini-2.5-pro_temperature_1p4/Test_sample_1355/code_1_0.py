import math

def double_factorial(n):
    """Computes the double factorial of an integer n."""
    if n < 0:
        return 0
    if n == 0 or n == 1:
        return 1
    result = 1
    for i in range(n, 1, -2):
        result *= i
    return result

def main():
    """
    Calculates the ratio between the fourth statistical moment of the dimensionless 
    conductance and its average value for a disordered Majorana wire at criticality.
    """
    print("This script calculates the ratio <g^4> / <g> for a disordered Majorana wire at the critical point.")
    print("The universal formula for the n-th moment is: <g^n> = (2n-1)!! / (2n)!!\n")

    # --- Calculate the average value <g> for n=1 ---
    n_avg = 1
    avg_num = double_factorial(2 * n_avg - 1)
    avg_den = double_factorial(2 * n_avg)
    print(f"First, we calculate the average conductance <g> (for n=1):")
    print(f"<g> = (2*1-1)!! / (2*1)!! = {avg_num}!! / {avg_den}!! = {avg_num} / {avg_den}\n")

    # --- Calculate the fourth moment <g^4> for n=4 ---
    n_4th = 4
    moment4_num = double_factorial(2 * n_4th - 1)
    moment4_den = double_factorial(2 * n_4th)
    print(f"Next, we calculate the fourth moment <g^4> (for n=4):")
    print(f"<g^4> = (2*4-1)!! / (2*4)!! = {2*n_4th-1}!! / {2*n_4th}!! = {moment4_num} / {moment4_den}\n")
    
    # --- Simplify the fraction for the fourth moment ---
    common_divisor = math.gcd(moment4_num, moment4_den)
    moment4_num_s = moment4_num // common_divisor
    moment4_den_s = moment4_den // common_divisor
    print(f"The fraction for <g^4> can be simplified from {moment4_num}/{moment4_den} to {moment4_num_s}/{moment4_den_s}.\n")

    # --- Calculate and print the final ratio ---
    # The ratio is (<g^4>) / (<g>) = (moment4_num / moment4_den) / (avg_num / avg_den)
    final_ratio_num = moment4_num_s * avg_den
    final_ratio_den = moment4_den_s * avg_num
    
    # Simplify the final ratio
    final_common_divisor = math.gcd(final_ratio_num, final_ratio_den)
    final_num = final_ratio_num // final_common_divisor
    final_den = final_ratio_den // final_common_divisor
    
    print("Finally, we compute the ratio <g^4> / <g>:")
    print(f"<g^4> / <g> = ({moment4_num_s}/{moment4_den_s}) / ({avg_num}/{avg_den}) = {final_num}/{final_den}")

if __name__ == "__main__":
    main()