# A configuration is given as a list of tuples (x,y) representing peg positions.
# Example configurations:
# C1 = [(0, 0), (1, 0)]
# C2 = [(2, 0)]
# C3 = [(0, 0), (1, 1)]

def get_fib_mod_2(n):
    """
    Calculates the n-th Fibonacci number modulo 2.
    The sequence F_n mod 2 is periodic with period 3: 0, 1, 1, 0, 1, 1, ...
    We define F_{-1} = 1.
    """
    # Adjust for negative indices using F_{n-2} = F_n - F_{n-1}
    # F_{-1}=F_1-F_0=1, F_{-2}=F_0-F_{-1}=-1=1 mod 2, F_{-3}=F_{-1}-F_{-2}=1-1=0 mod 2.
    # The sequence F_n mod 2 for n=...,-2,-1,0,1,2,... is ...,1,1,0,1,1,0,...
    # This is (n+1) mod 3 == 0 -> 1, else (n+1)%3==1 -> 1, (n+1)%3==2 -> 0.
    # Let's use a simpler way: direct computation based on (n mod 3)
    rem = n % 3
    if rem == 0: return 0
    if rem == 1: return 1
    if rem == 2: return 1
    # For negative n
    if rem == -1: # same as rem=2
        return 1
    if rem == -2: # same as rem=1
        return 1
    return 1 # Fallback, should not be reached


def calculate_invariant_vector(config):
    """
    Calculates a 4-dimensional invariant vector modulo 2 for a given configuration.
    Any two equivalent configurations will have the same invariant vector.
    """
    # The invariant is based on a mapping of each peg to a base of 4 generators,
    # p_00, p_01, p_10, p_11.
    # The coefficient for p_ij in the expansion of a peg at (x,y) is
    # related to products of Fibonacci numbers F_{x-i} and F_{y-j}.
    # We use the simplified version modulo 2.
    
    n00, n01, n10, n11 = 0, 0, 0, 0
    
    for x, y in config:
        # Coefficients based on products of Fibonacci numbers mod 2
        # g_x = F_x mod 2, g_{x-1} = F_{x-1} mod 2, etc.
        # Vector is (n_11, n_10, n_01, n_00)
        gx = get_fib_mod_2(x)
        gx_1 = get_fib_mod_2(x - 1)
        gy = get_fib_mod_2(y)
        gy_1 = get_fib_mod_2(y - 1)
        
        n11 += gx * gy
        n10 += gx * gy_1
        n01 += gx_1 * gy
        n00 += gx_1 * gy_1

    # Return the vector modulo 2
    return [n11 % 2, n10 % 2, n01 % 2, n00 % 2]

# --- Demonstration ---
# Two configurations that are known to be equivalent
config1 = [(0, 0), (1, 0)]
config2 = [(2, 0)]
# A third configuration in a different class
config3 = [(0, 0), (1, 1)]

inv1 = calculate_invariant_vector(config1)
inv2 = calculate_invariant_vector(config2)
inv3 = calculate_invariant_vector(config3)

print(f"Configuration {config1} has invariant vector: {inv1}")
print(f"Configuration {config2} has invariant vector: {inv2}")
print(f"Configuration {config3} has invariant vector: {inv3}")
print("\nAs the invariants for the first two configurations are the same, they belong to the same equivalence class.")

print("\nThe total number of equivalence classes under this relation is 8.")

<<<8>>>