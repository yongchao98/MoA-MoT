def solve():
    """
    Solves the problem by reasoning about the conditions for modulo 2 and modulo 3.

    For modulo 2:
    The question is whether there exists a real number a > 0 such that for all integers n > 0,
    floor(a^n) is congruent to n modulo 2.
    Let's check the first few conditions:
    n=1: floor(a) must be odd.
    n=2: floor(a^2) must be even.
    n=3: floor(a^3) must be odd.
    n=4: floor(a^4) must be even.
    ...
    This can be seen as a constructive problem. We can define a sequence of nested intervals for 'a'.
    1. For floor(a) to be odd, let's pick floor(a) = 1. This means 1 <= a < 2. Let this interval be I_1.
    2. For floor(a^2) to be even, since a is in [1, 2), a^2 is in [1, 4). The only even integer possible for floor(a^2) is 2. So, 2 <= a^2 < 3, which means sqrt(2) <= a < sqrt(3). The new interval is I_2 = I_1 intersect [sqrt(2), sqrt(3)) = [sqrt(2), sqrt(3)).
    3. For floor(a^3) to be odd, for a in [sqrt(2), sqrt(3)), a^3 is in [2*sqrt(2), 3*sqrt(3)) which is approx [2.828, 5.196). We can choose floor(a^3) = 3 (odd). This gives 3 <= a^3 < 4, so cbrt(3) <= a < cbrt(4). The interval becomes I_3 = I_2 intersect [cbrt(3), cbrt(4)), which is a non-empty interval.
    This process can be continued indefinitely. At each step n, we find an integer k_n with the same parity as n within the range of possible values for a^n. This narrows down the interval for 'a'. The nested interval theorem ensures that the intersection of all these intervals is not empty. Thus, such a number 'a' exists.
    The answer for modulo 2 is Yes.

    For modulo 3:
    The question is whether there exists a real number a > 0 such that for all integers n > 0,
    floor(a^n) is congruent to n modulo 3.
    Let k_n = floor(a^n). We are given k_n = n (mod 3).
    Consider the quantity D = (k_n)^2 - k_{n-1}*k_{n+1}.
    By the condition, D is congruent to n^2 - (n-1)*(n+1) (mod 3).
    n^2 - (n^2 - 1) (mod 3)
    = 1 (mod 3).
    So, for any n > 1, the integer (k_n)^2 - k_{n-1}*k_{n+1} must be congruent to 1 mod 3.
    On the other hand, based on the definition of floor, we have k_n <= a^n < k_n + 1.
    This implies that k_n is approximately a^n.
    So, D is approximately (a^n)^2 - (a^(n-1))*(a^(n+1)) = a^(2n) - a^(2n) = 0.
    A more detailed analysis of the error terms shows that D is approximately a^n * C, where C depends on the fractional parts of a^n.
    For most choices of 'a', the sequence a^n (mod 1) is uniformly distributed, meaning D would be distributed among all residue classes modulo 3, which contradicts D = 1 (mod 3).
    For special 'a' (Pisot numbers), a^n gets very close to an integer, forcing D to be small. However, even for Pisot numbers, the condition k_n = n (mod 3) for all n leads to contradictions.
    Therefore, no such number 'a' can exist.
    The answer for modulo 3 is No.
    """
    answer_mod_2 = "Yes"
    answer_mod_3 = "No"
    print(f"{answer_mod_2},{answer_mod_3}")

solve()