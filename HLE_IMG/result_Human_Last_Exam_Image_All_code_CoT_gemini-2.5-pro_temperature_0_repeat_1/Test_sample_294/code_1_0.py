import math

def solve():
    """
    This function finds the minimum positive integer N such that S(N) is not empty.
    S(N) is the set of pairs of positive integers (n1, n2) with n1 <= N and n2 <= N
    for which the group GG_{n1,n2}(r) is infinite for some r > 0.

    The group is infinite if and only if:
    1/n1 + 1/n2 + 1/lcm(n1, n2) <= 1.

    The algorithm iterates through N = 1, 2, 3, ... and for each N, checks if there
    exists a pair (n1, n2) with 1 <= n1, n2 <= N that satisfies the inequality.
    The first N for which such a pair is found is the answer.
    """

    def lcm(a, b):
        """Computes the least common multiple of a and b."""
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

    N = 1
    while True:
        # For a given N, search for a pair (n1, n2) where 1 <= n1, n2 <= N
        # that satisfies the condition for an infinite group.
        found_pair = None
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                l = lcm(n1, n2)
                # The condition for an infinite group is 1/n1 + 1/n2 + 1/l <= 1.
                # We use integer arithmetic to avoid floating point inaccuracies:
                # (n2*l + n1*l + n1*n2) <= (n1 * n2 * l)
                if (n2 * l + n1 * l + n1 * n2) <= (n1 * n2 * l):
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break

        if found_pair:
            n1, n2 = found_pair
            l = lcm(n1, n2)
            val = 1/n1 + 1/n2 + 1/l
            
            print("To find the minimum N, we test values of N starting from 1.")
            print("For N=1 and N=2, no pair (n1, n2) satisfies the condition 1/n1 + 1/n2 + 1/lcm(n1,n2) <= 1. So, S(1) and S(2) are empty.")
            print(f"\nNow, we test N={N}.")
            print(f"We found a pair (n1, n2) = ({n1}, {n2}) where both numbers are less than or equal to {N}.")
            print("The condition for the group to be infinite is 1/n1 + 1/n2 + 1/lcm(n1, n2) <= 1.")
            print(f"For the pair ({n1}, {n2}), the least common multiple is lcm({n1}, {n2}) = {l}.")
            print("We check the sum:")
            print(f"1/{n1} + 1/{n2} + 1/{l} = {1/n1:.6f} + {1/n2:.6f} + {1/l:.6f} = {val}")
            print(f"Since {val} <= 1, the condition is met.")
            print(f"This means the pair ({n1}, {n2}) is in S({N}), so S({N}) is not empty.")
            print(f"Therefore, the minimum N for which S(N) is not empty is {N}.")
            
            print(f"\n<<<{N}>>>")
            return

        N += 1

if __name__ == '__main__':
    solve()