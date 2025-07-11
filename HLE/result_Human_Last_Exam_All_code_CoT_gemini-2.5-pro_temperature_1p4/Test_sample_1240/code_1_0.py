import math

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def solve_2a_mod_d(d):
    """Finds solutions for 2*a = 0 (mod d) in {0, 1, ..., d-1}."""
    solutions = []
    for a in range(d):
        if (2 * a) % d == 0:
            solutions.append(a)
    return solutions

def analyze_question_a():
    """
    Q: Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?
    """
    print("--- Analyzing Question (a) ---")
    print("Is it possible for R_2(M) of a d-neighbor of Z^12 to be type A_11?")
    n = 12
    
    print("An A_11 root system consists of vectors e_i - e_j for 1 <= i,j <= 12.")
    print("For a vector v to be in M, we need v . w = 0 (mod d).")
    
    print("\nStep 1: Find conditions on w for e_i - e_j to be in M.")
    print("  (e_i - e_j) . w = w_i - w_j. So we need w_i - w_j = 0 (mod d) for all i, j.")
    print("  This means all components of w must be congruent: w_1 = w_2 = ... = w_12 = c (mod d).")
    
    print("\nStep 2: Find conditions for e_i + e_j to be excluded from M.")
    print("  (e_i + e_j) . w = w_i + w_j. So we need w_i + w_j != 0 (mod d).")
    print("  With w_i = c, this becomes c + c = 2*c != 0 (mod d).")

    print("\nStep 3: Check the index condition.")
    print("  The sublattice M has index d in Z^n only if gcd(w_1, ..., w_n, d) = 1.")
    print("  With w_i = c, this is gcd(c, ..., c, d) = gcd(c, d) = 1.")
    
    print("\nStep 4: Find a valid c and d.")
    print("  Let's choose the simplest c satisfying gcd(c, d) = 1, which is c = 1.")
    c = 1
    print(f"  With c = {c}, the condition becomes 2 * {c} != 0 (mod d), which is 2 != 0 (mod d).")
    print("  This means d cannot be 1 or 2. Let's pick d = 3.")
    d = 3
    
    print(f"\nConstruction: n=12, d=3, w = ({c}, ..., {c}).")
    print(f"  Check e_i - e_j: w_i - w_j = {c} - {c} = 0.  0 mod 3 = 0. This works.")
    print(f"  Check e_i + e_j: w_i + w_j = {c} + {c} = 2.  2 mod 3 != 0. This works.")
    print(f"  Check index: gcd({c}, {d}) = {math.gcd(c, d)}. This is 1, so it works.")

    print("\nConclusion for (a): A valid construction exists. The answer is Yes.")
    return "Yes"

def analyze_question_b():
    """
    Q: Can R_2(M) of a d-neighbor of Z^15 contain a D_7 component?
    """
    print("\n--- Analyzing Question (b) ---")
    print("Can R_2(M) of a d-neighbor of Z^15 contain a D_7 component?")
    n = 15
    k = 7 # D_7
    
    print("A D_k component on indices S requires e_i +/- e_j to be in M for all i,j in S.")
    print("  (e_i - e_j) . w = 0 (mod d) => w_i = w_j = a (mod d) for i,j in S.")
    print("  (e_i + e_j) . w = 0 (mod d) => w_i + w_j = 2a = 0 (mod d).")
    print("Let's try to construct such a w and d.")

    print("\nStep 1: Choose d and find a.")
    print("  Let's try d = 2. The condition is 2*a = 0 (mod 2).")
    a_test = 1
    print(f"  For d=2, {2*a_test} = 2, and 2 mod 2 = 0. So any 'a' works. Let's pick a = 1.")
    d = 2
    a = 1
    
    print("\nStep 2: Construct w to isolate the component.")
    print("  Let S = {1, ..., 7}. We set w_i = a = 1 for i in S.")
    print("  To make D_7 a separate component, roots connecting S and {8, ..., 15} must not be in M.")
    print("  Let i be in S and j not in S. We need w_i +/- w_j != 0 (mod 2).")
    print(f"  With w_i = {a}, we need {a} +/- w_j != 0 (mod 2). So 1 +/- w_j must be odd.")
    print("  This forces w_j to be even, i.e., w_j = 0 (mod 2).")
    
    w = [a] * k + [0] * (n - k)
    print(f"\nProposed construction: d=2, w = {w}")

    print("\nStep 3: Verify the construction.")
    w_i = 1 # in S
    w_j = 0 # not in S
    print(f"  Isolation check: w_i + w_j = {w_i}+{w_j} = 1 (mod 2). Correct.")
    print(f"                   w_i - w_j = {w_i}-{w_j} = 1 (mod 2). Correct.")
    gcd_val = gcd_list(w + [d])
    print(f"  Index condition: gcd(of w components and d) = {gcd_val}. Correct.")

    print("\nThis construction yields a root system with a D_7 component. In fact, it is D_7 + D_8.")
    print("Conclusion for (b): A valid construction exists. The answer is yes.")
    return "yes"

def analyze_question_c():
    """
    Q: For n = 18, d = 5, is it possible for R_2(M) to have >1 D-component?
    """
    print("\n--- Analyzing Question (c) ---")
    print("For n=18, d=5, can R_2(M) have more than one D-type component?")
    d = 5

    print("As shown before, a D-type component on indices S requires w_i = a (mod d) for i in S,")
    print(f"and also that 2*a = 0 (mod d). Let's solve this for d = {d}.")
    
    solutions = solve_2a_mod_d(d)
    a_sol = solutions[0]
    print(f"  The equation 2*a = 0 (mod 5) has only one solution in Z_5: a = {a_sol}.")
    print("  This is because gcd(2, 5) = 1, so we can divide by 2.")

    print("\nThis result implies that for any D-type component to exist, the corresponding")
    print("glue vector components w_i must all be 0 (mod 5).")

    print("\nNow, assume R_2(M) has two D-components, D_k on indices S_k and D_l on S_l.")
    print("  This means w_i = 0 (mod 5) for all i in S_k.")
    print("  And w_j = 0 (mod 5) for all j in S_l.")

    print("\nLet's check a root `v = e_i + e_j` connecting the two components (i in S_k, j in S_l).")
    w_i, w_j = 0, 0
    dot_product = w_i + w_j
    print(f"  v . w = w_i + w_j = {w_i} + {w_j} = {dot_product}.")
    print(f"  The result modulo {d} is {dot_product % d}.")
    print("  Since the result is 0, this connecting root is in M.")
    print("  This contradicts the assumption that D_k and D_l are separate components.")
    print("  They are merged into a single, larger D_{k+l} component.")
    
    print("\nConclusion for (c): It is not possible to have more than one D-component. The answer is no.")
    return "no"

def solve():
    """Main function to solve the user's question."""
    ans_a = analyze_question_a()
    ans_b = analyze_question_b()
    ans_c = analyze_question_c()
    
    final_answer = f"(a) [{ans_a}]; (b) [{ans_b}]; (c) [{ans_c}]."
    print("\n" + "="*50)
    print("Final Summary:")
    print(final_answer)
    print(f"\n<<<(a) [Yes]; (b) [yes]; (c) [no]>>>")

solve()