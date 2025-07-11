#
# This program determines how many integers in a given range can be expressed
# in the form n = x^3 + 2y^3 + 4z^3 - 6xyz.
#
# The solution is based on a property of the expression modulo 9.
#
# 1. A number n can be expressed in the given form only if n mod 9 is 0, 1, or 8.
#    This is a necessary condition. For the large numbers in the problem's range,
#    this condition is also sufficient.
#
# 2. The problem reduces to counting the integers n in the range
#    [10^18, 10^18 + 10000] that satisfy this modular condition.
#
# 3. We determine the modular value of the start of the range: 10^18 mod 9.
#    Since 10 = 1 (mod 9), 10^18 = 1^18 = 1 (mod 9).
#
# 4. We then count the number of integers k in [0, 10000] such that (1 + k) mod 9
#    is 0, 1, or 8. This is equivalent to k mod 9 being 8, 0, or 7.
#
def solve():
    """
    Counts the number of integers k in [0, 10000] where k mod 9 is 0, 7, or 8.
    """
    k_min = 0
    k_max = 10000

    # Count numbers where k === 0 (mod 9)
    # k = 9*m; 0 <= 9*m <= 10000  =>  0 <= m <= 1111.11...
    # m can be 0, 1, ..., 1111.
    count_0_mod_9 = 1111 - 0 + 1

    # Count numbers where k === 7 (mod 9)
    # k = 9*m + 7; 0 <= 9*m + 7 <= 10000  =>  -7/9 <= m <= 9993/9 = 1110.33...
    # m can be 0, 1, ..., 1110.
    count_7_mod_9 = 1110 - 0 + 1

    # Count numbers where k === 8 (mod 9)
    # k = 9*m + 8; 0 <= 9*m + 8 <= 10000  =>  -8/9 <= m <= 9992/9 = 1110.22...
    # m can be 0, 1, ..., 1110.
    count_8_mod_9 = 1110 - 0 + 1

    total_count = count_0_mod_9 + count_7_mod_9 + count_8_mod_9

    print("The number of integers is found by summing the counts for each valid modular case:")
    print(f"Count for k = 0 (mod 9): {count_0_mod_9}")
    print(f"Count for k = 7 (mod 9): {count_7_mod_9}")
    print(f"Count for k = 8 (mod 9): {count_8_mod_9}")
    print(f"\nTotal count = {count_0_mod_9} + {count_7_mod_9} + {count_8_mod_9} = {total_count}")

solve()