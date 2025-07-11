import math

def final_calculation():
    """
    This function calculates the cardinality of [F * G](.) for n=8.
    It follows the plan of using Burnside's Lemma on the set A x Hom(A,A).
    """

    n = 8
    A_nums = list(range(1, n + 1))
    # 'inf' represents the infinity element of the monoid A
    A_vals = A_nums + ['inf']
    A_size = len(A_vals)

    # Memoization cache for gcd to speed up calculations
    _gcd_cache = {}
    def gcd(a, b):
        # Canonical order for cache key
        if not isinstance(a, str) and not isinstance(b, str) and a > b:
            a,b = b,a
        elif isinstance(a, str) and not isinstance(b, str): # inf should be first
            a, b = b, a
            
        if (a, b) in _gcd_cache:
            return _gcd_cache[(a,b)]
        
        if a == 'inf': res = b
        elif b == 'inf': res = a
        else: res = math.gcd(a,b)
        
        _gcd_cache[(a,b)] = res
        return res

    # --- Step 1: Find all homomorphisms f: A -> A ---
    homomorphisms = []
    
    # Recursive function to find homomorphisms by building them element by element
    def find_homs_recursive(k, f):
        if k > n:
            f_final = f.copy()
            f_final['inf'] = 'inf'
            homomorphisms.append(f_final)
            return

        for fk_candidate in A_vals:
            f[k] = fk_candidate
            possible = True
            # Check if the new f(k) is consistent with the homomorphism property for previously set values
            for i in range(1, k):
                if f[gcd(i, k)] != gcd(f[i], fk_candidate):
                    possible = False
                    break
            if possible:
                find_homs_recursive(k + 1, f)
        
        # Backtrack
        if k in f:
            del f[k]

    # Initialize the recursion: f(1) can be any value in A
    f_base = {}
    for f1 in A_vals:
        f_base[1] = f1
        find_homs_recursive(2, f_base)
        
    Y = homomorphisms
    Y_size = len(Y)
    
    total_fixed_points = 0
    
    print(f"Plan: Calculate Sum(|X^a| * |Y^a|) for a in A, then divide by |A|={A_size}.")
    print(f"Found |Y| = {Y_size} homomorphisms.")
    print("-" * 20)
    print("a\t|X^a|\t|Y^a|\tTerm")
    print("-" * 20)
    
    # --- Step 2 & 3: Loop through a, calculate terms, and sum ---
    all_terms = []
    for a in A_vals:
        # Calculate |X^a| = number of divisors of a in A
        X_a_size = 0
        for x in A_vals:
            if gcd(a, x) == x:
                X_a_size += 1
                
        # Calculate |Y^a| = number of fixed homomorphisms
        Y_a_size = 0
        for f in Y:
            is_fixed = True
            for b in A_vals:
                if f[b] != f[gcd(a, b)]:
                    is_fixed = False
                    break
            if is_fixed:
                Y_a_size += 1
        
        term = X_a_size * Y_a_size
        total_fixed_points += term
        all_terms.append(term)
        
        a_str = str(a) if a != 'inf' else 'inf'
        print(f"{a_str:<4}\t{X_a_size:<5}\t{Y_a_size:<5}\t{term}")

    print("-" * 20)
    
    # --- Step 4: Final calculation and output ---
    equation_str = " + ".join(map(str, all_terms))
    print(f"Sum of terms: {equation_str} = {total_fixed_points}")
    
    result = total_fixed_points // A_size
    
    print(f"Final cardinality = ({total_fixed_points}) / {A_size} = {result}")

final_calculation()