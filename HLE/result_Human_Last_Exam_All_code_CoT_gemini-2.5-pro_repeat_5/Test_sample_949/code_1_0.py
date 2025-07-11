import math

def calculate_sums(k):
    """
    Calculates the LHS and RHS sums for the 6 statements
    for the string Sk = (()()...()) with k inner pairs.
    """
    
    # For Sk, there is one outer pair (x0) and k inner pairs (xi).
    # L(x0) = 2k + 2, D(x0) = 2
    # L(xi) = 2, D(xi) = 1
    
    results = []

    # Statement 1: sum log L = O(sum log D)
    lhs1 = math.log(2 * k + 2) + k * math.log(2)
    rhs1 = math.log(2) + k * math.log(1) if k * math.log(1) == 0 else math.log(2) # log(1)=0
    results.append(lhs1 / rhs1 if rhs1 != 0 else float('inf'))

    # Statement 2: sum loglog L = O(sum loglog D)
    # Sum over well-defined terms. log(log(v)) is defined for v > 1.
    # Base is assumed to be math.e > 2, so log(2)<1, log(log(2)) is undefined.
    lhs2 = 0
    if 2 * k + 2 > math.e: # log(v)>1
      lhs2 += math.log(math.log(2*k+2))
    # No L(xi) or D(x) term is > e, so only one term on LHS is summed
    rhs2 = 0 # No D(x) is > e
    results.append(lhs2 / rhs2 if rhs2 != 0 else float('inf'))


    # Statement 3: sum log^5 L = O(sum log^5 D)
    lhs3 = math.pow(math.log(2 * k + 2), 5) + k * math.pow(math.log(2), 5)
    rhs3 = math.pow(math.log(2), 5) # log(1)=0
    results.append(lhs3 / rhs3 if rhs3 != 0 else float('inf'))

    # Statement 4: sum 2^sqrt(log L) = O(sum 2^sqrt(log D))
    # We check if LHS <= C * RHS for some C. Here we just show the ratio.
    # The O() is inside, so we assume C=2 for demonstration.
    lhs4 = 2**math.sqrt(math.log(2 * k + 2)) + k * 2**math.sqrt(math.log(2))
    rhs4 = 2**(2*math.sqrt(math.log(2))) + k * 2**(2*math.sqrt(math.log(1)))
    results.append(lhs4 / rhs4 if rhs4 != 0 else float('inf'))

    # Statement 5: sum L^0.1 = O(sum D^0.11)
    lhs5 = (2 * k + 2)**0.1 + k * (2**0.1)
    rhs5 = (2**0.11) + k * (1**0.11)
    results.append(lhs5 / rhs5 if rhs5 != 0 else float('inf'))
    
    # Statement 6: sum L^0.25 = O(sum D^0.5)
    lhs6 = (2 * k + 2)**0.25 + k * (2**0.25)
    rhs6 = (2**0.5) + k * (1**0.5)
    results.append(lhs6 / rhs6 if rhs6 != 0 else float('inf'))

    return results

def main():
    """
    Main function to run the demonstration.
    """
    print("This script analyzes the 6 statements for the string family Sk = (()()...()).")
    print("It calculates the ratio LHS_sum / RHS_sum for growing k.")
    print("For FALSE statements, this ratio should grow indefinitely.")
    print("For TRUE statements, this ratio should be bounded.\n")
    
    print(" k   |  S1 Ratio |  S2 Ratio |  S3 Ratio |  S4 Ratio |  S5 Ratio |  S6 Ratio ")
    print("--------------------------------------------------------------------------------")
    
    for k_power in range(1, 6):
        k = 10**k_power
        ratios = calculate_sums(k)
        print(f"{k:<6d} | {ratios[0]:>10.2f} | {ratios[1]:>10.2f} | {ratios[2]:>10.2f} | {ratios[3]:>10.2f} | {ratios[4]:>10.2f} | {ratios[5]:>10.2f}")

if __name__ == "__main__":
    main()
