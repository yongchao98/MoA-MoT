import fractions

def find_min_N():
    """
    Finds the minimum positive integer N such that there exists a pair (n1, n2)
    with 1 <= n1, n2 <= N for which the group GG_{n1, n2}(r) can be infinite.
    
    The condition for the group to be infinite is 1/n1 + 1/n2 + 1/k <= 1,
    where k is the order of the composition of the two rotations, which is
    the denominator of the fraction 1/n1 + 1/n2.
    """
    N = 1
    while True:
        print(f"Testing N = {N}...")
        found_pair = None
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # Trivial cases: if n1 or n2 is 1, the group is finite (C_n).
                if n1 == 1 or n2 == 1:
                    continue

                # Calculate k
                try:
                    frac_sum = fractions.Fraction(1, n1) + fractions.Fraction(1, n2)
                    k = frac_sum.denominator
                except ZeroDivisionError:
                    continue
                
                # Check the condition 1/n1 + 1/n2 + 1/k <= 1
                # We use integer arithmetic to avoid floating point inaccuracies.
                # Equivalent to n2*k + n1*k + n1*n2 <= n1*n2*k
                if (n2 * k + n1 * k + n1 * n2) <= (n1 * n2 * k):
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            frac_sum = fractions.Fraction(1, n1) + fractions.Fraction(1, n2)
            k = frac_sum.denominator
            total_sum_frac = frac_sum + fractions.Fraction(1, k)
            
            print(f"For N = {N}, found a working pair (n1, n2) = ({n1}, {n2}).")
            print("The condition for the group to be infinite is 1/n₁ + 1/n₂ + 1/k ≤ 1.")
            print(f"For this pair, k (the order of the composite rotation) is {k}.")
            print(f"Checking the condition for ({n1}, {n2}):")
            print(f"1/{n1} + 1/{n2} + 1/{k} = {frac_sum.numerator}/{frac_sum.denominator} + 1/{k} = {total_sum_frac}")
            print(f"Since {total_sum_frac} <= 1, the condition is met.")
            print(f"The minimum value of N is {N}.")
            return N
        else:
            print(f"For N = {N}, no such pair was found. S({N}) is empty.")
            N += 1

if __name__ == '__main__':
    min_N = find_min_N()
    # The final answer in the specified format
    # print(f"<<<{min_N}>>>")