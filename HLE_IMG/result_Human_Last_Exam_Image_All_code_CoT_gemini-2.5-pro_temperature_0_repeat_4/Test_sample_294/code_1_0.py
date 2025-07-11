import fractions

def find_min_N():
    """
    Finds the minimum positive integer N such that there exists a pair of positive integers (n1, n2)
    with n1 <= N and n2 <= N for which the group GG_{n1,n2}(r) is infinite for some r > 0.

    The condition for the group to be infinite is 1/n1 + 1/n2 < 1.
    """
    N = 1
    while True:
        # For a given N, search for a pair (n1, n2) satisfying the conditions.
        # We must have n1 >= 2 and n2 >= 2, because if either is 1, the group is finite.
        found_pair = False
        solution_pair = None

        # The function 1/n1 + 1/n2 decreases as n1 or n2 increase.
        # We are looking for the smallest N, so we check pairs with small n1, n2 first.
        if N >= 2:
            for n1 in range(2, N + 1):
                for n2 in range(2, N + 1):
                    # Use fractions for exact arithmetic
                    sum_of_inverses = fractions.Fraction(1, n1) + fractions.Fraction(1, n2)
                    if sum_of_inverses < 1:
                        found_pair = True
                        solution_pair = (n1, n2)
                        break
                if found_pair:
                    break
        
        if found_pair:
            n1, n2 = solution_pair
            sum_val = fractions.Fraction(1, n1) + fractions.Fraction(1, n2)
            
            print("The condition for the group GG_{n₁,n₂}(r) to be infinite for some r > 0 is 1/n₁ + 1/n₂ < 1.")
            print("We seek the smallest integer N such that a pair (n₁, n₂) exists with n₁ ≤ N, n₂ ≤ N satisfying this.")
            print("\nTesting N=2:")
            print("The only pair with n₁, n₂ ≥ 2 is (2, 2).")
            print("For (2, 2), the sum is 1/2 + 1/2 = 1. This is not less than 1. So S(2) is empty.")
            
            print(f"\nTesting N={N}:")
            print(f"A pair satisfying the condition is (n₁, n₂) = ({n1}, {n2}).")
            print("Let's check the final equation for this pair:")
            print(f"1/{n1} + 1/{n2} = {sum_val.numerator}/{sum_val.denominator}")
            print(f"Since {sum_val.numerator}/{sum_val.denominator} < 1, the condition is met for the pair ({n1}, {n2}).")
            print(f"For this pair, we need N ≥ n₁ and N ≥ n₂. So, N ≥ {n1} and N ≥ {n2}, which means N ≥ {max(n1, n2)}.")
            print(f"Since no solution exists for N < {N}, the minimum value for N is {N}.")
            
            return N
            
        N += 1

min_N = find_min_N()
# The final answer is the value of min_N
# print(f"\n<<< {min_N} >>>")