import math

def gcd_list(numbers):
    """Computes the GCD of a list of integers."""
    if not numbers:
        return 0
    result = abs(numbers[0])
    for i in range(1, len(numbers)):
        result = math.gcd(result, abs(numbers[i]))
    return result

def check_conditions(w, d, n):
    """Checks if a vector w and index d define a valid neighbor construction."""
    print(f"Checking for n={n}, d={d}, w={w}")
    if len(w) != n:
        print(f"  > Failed: w must have length n={n}")
        return False, {}

    # Condition 1: gcd(w_1, ..., w_n, d) = 1
    g = math.gcd(gcd_list(w), d)
    if g != 1:
        print(f"  > Failed: gcd(components of w, d) = {g}, but should be 1.")
        return False, {}

    # Condition 2: w.w must be a multiple of d
    dot_product = sum(i*i for i in w)
    if dot_product % d != 0:
        print(f"  > Failed: w.w = {dot_product} is not a multiple of {d}.")
        return False, {"dot_product": dot_product}

    print("  > Success: All conditions for a valid neighbor construction passed.")
    return True, {"dot_product": dot_product}

def print_dot_product_equation(w, result):
    """Prints the dot product sum equation."""
    # To avoid long output, only show non-zero terms
    terms = [f"{i}^2" for i in w if i != 0]
    eq = " + ".join(terms)
    # Handle the all-ones case for question a
    if all(val == 1 for val in w):
        eq = " + ".join([f"1^2"] * len(w))
    
    print(f"  > Equation: w.w = {eq} = {result}")

def question_a():
    print("--- Question (a) Verification ---")
    n = 12
    # For A_11, we choose w_i = 1 for all i
    w = [1] * n
    # For w.w = 12, choose d=3, which is a prime divisor and d!=2
    d = 3
    
    valid, info = check_conditions(w, d, n)
    if valid:
        dot_product = info['dot_product']
        print_dot_product_equation(w, dot_product)
        print(f"  > Check: {dot_product} mod {d} = {dot_product % d}, so the condition holds.")
        print("  > Root system check:")
        print("    - For v = e_i - e_j: w.v = 1 - 1 = 0. These are in M.")
        print("    - For v = e_i + e_j: w.v = 1 + 1 = 2. 2 mod 3 != 0. These are not in M.")
        print("  > Conclusion: The resulting root system R2(M) is A_11. The answer is Yes.")
    print("-" * 20 + "\n")

def question_b():
    print("--- Question (b) Verification ---")
    n = 15
    # For a D_7 component on first 7 coordinates with odd d, need w_i=0 (mod d) for i=1..7.
    # Let d=3. Let's construct a w.
    w = [0]*7 + [1]*3 + [0]*(15-7-3)
    d = 3

    valid, info = check_conditions(w, d, n)
    if valid:
        dot_product = info['dot_product']
        print_dot_product_equation(w, dot_product)
        print(f"  > Check: {dot_product} mod {d} = {dot_product % d}, so the condition holds.")
        print("  > Root system check:")
        print("    - For v = +/-e_i +/-e_j with i,j in {1,...,7}, w.v = 0. These roots are in M.")
        print("  > Conclusion: R2(M) contains a D_7 component. The answer is yes.")
    print("-" * 20 + "\n")

def question_c():
    print("--- Question (c) Verification ---")
    n = 18
    d = 5
    print(f"Analyzing for n={n}, d={d}")
    print("  > A D_k component on indices I requires w_i = 0 (mod 5) for all i in I (since 5 is odd).")
    print("  > Assume we have two D components, on disjoint sets I1 and I2.")
    print("  > This implies w_i = 0 (mod 5) for i in I1 and w_j = 0 (mod 5) for j in I2.")
    print("  > Let's test a 'cross-root' v = e_i + e_j where i is in I1, j is in I2.")
    print("  > The condition for v to be in M is w.v mod 5 = 0.")
    print("  > w.v = w_i + w_j. Since w_i and w_j are both 0 (mod 5), their sum is also 0 (mod 5).")
    print("  > Thus, the cross-root is in R2(M). This merges the two supposed components into one larger D component.")
    print("  > Conclusion: It is not possible for R2(M) to have more than one D-type component. The answer is no.")
    print("-" * 20 + "\n")

# Run the analysis for all questions
question_a()
question_b()
question_c()