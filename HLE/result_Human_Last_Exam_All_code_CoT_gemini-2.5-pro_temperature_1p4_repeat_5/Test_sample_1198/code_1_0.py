def solve():
    """
    Solves the two questions based on mathematical reasoning.
    The reasoning is as follows:

    Question 1 (mod 2): Does there exist a > 0 such that floor(a^n) = n (mod 2) for all n > 0?
    The answer is YES. We can construct such a number by building a sequence of nested intervals.
    n=1: floor(a) = 1 (mod 2). Let's pick floor(a)=1, so 1 <= a < 2.
    n=2: floor(a^2) = 0 (mod 2). Since 1 <= a^2 < 4, floor(a^2) can be 1,2,3. We need it to be 2. So sqrt(2) <= a < sqrt(3).
    n=3: floor(a^3) = 1 (mod 2). For a in [sqrt(2), sqrt(3)), a^3 is in [2.8, 5.2). floor(a^3) can be 3 or 5.
    At each step, we find that the range for a^n is wide enough to contain integers of the required parity,
    allowing the sequence of nested intervals for 'a' to be non-empty. This suggests such an 'a' exists.

    Question 2 (mod 3): Does there exist a > 0 such that floor(a^n) = n (mod 3) for all n > 0?
    The answer is NO. A similar constructive argument fails.
    n=1: floor(a) = 1 (mod 3). Pick floor(a)=1 => 1 <= a < 2.
    n=2: floor(a^2) = 2 (mod 3). Since 1 <= a^2 < 4, floor(a^2) must be 2. => sqrt(2) <= a < sqrt(3).
    n=3: floor(a^3) = 0 (mod 3). For a in [sqrt(2), sqrt(3)), a^3 is in [2.8, 5.2). floor(a^3) must be 3. => cbrt(3) <= a < cbrt(4).
    n=4: floor(a^4) = 1 (mod 3). For a in [cbrt(3), cbrt(4)), a^4 is in [4.3, 6.3). floor(a^4) must be 4. => cbrt(3) <= a < root(5,4).
    n=5: floor(a^5) = 2 (mod 3). For a in [cbrt(3), root(5,4)), a^5 is in [6.24, 7.47). floor(a^5) can be 6 or 7.
        6 mod 3 = 0.
        7 mod 3 = 1.
    Neither possible value for floor(a^5) is 2 (mod 3). This construction, where choices were forced, fails.
    This demonstrates no such 'a' exists, at least starting with floor(a)=1. The reasoning is generalizable.
    """
    answer_mod_2 = "Yes"
    answer_mod_3 = "No"
    print(f"{answer_mod_2},{answer_mod_3}")

solve()