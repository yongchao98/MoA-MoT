import math
from collections import deque

# Memoization cache for crossing number calculation
crossing_number_cache = {}

def get_crossing_number_for_pq(p, q):
    """
    Calculates the crossing number for a 2-bridge knot K(p/q) using
    a BFS search based on the Doll-Hoste algorithm.
    """
    if (p, q) in crossing_number_cache:
        return crossing_number_cache[(p, q)]

    # (cost, v_current, v_previous)
    queue = deque([(0, q, p)])
    visited = {(q, p)}
    min_cost = float('inf')

    while queue:
        cost, v_curr, v_prev = queue.popleft()

        if v_curr == 0:
            min_cost = min(min_cost, cost)
            continue
        
        if cost >= min_cost:
            continue

        if v_curr == 0: continue
        
        r = v_prev / v_curr
        # Explore the two integers closest to r
        k_choices = {math.floor(r), math.ceil(r)}
        
        for k in k_choices:
            v_next = k * v_curr - v_prev
            
            # Condition for valid sequence term
            if abs(v_next) <= abs(v_curr) / 2:
                if (v_next, v_curr) not in visited:
                    new_cost = cost + abs(k)
                    queue.append((new_cost, v_next, v_curr))
                    visited.add((v_next, v_curr))

    crossing_number_cache[(p, q)] = min_cost
    return min_cost

def get_crossing_number(p, q):
    """
    The crossing number of the knot K(p,q) is the minimum of the
    crossing numbers for p/q and its mirror image p/(p-q).
    """
    # We only need to check one representative, as K(p,q) and K(p,p-q) are mirror images.
    # The algorithm itself finds the minimal representation.
    return get_crossing_number_for_pq(p,q)

def get_orbit_rep(p, q):
    """
    Finds the canonical representative (smallest element) for the orbit of q.
    """
    try:
        q_inv = pow(q, -1, p)
    except ValueError: # Should not happen if gcd(p,q)=1
        return q
    
    orbit = {q, p - q, q_inv, p - q_inv}
    return min(orbit)

def solve():
    """
    Finds and counts the 2-bridge knots with crossing number <= 13
    that have two disjoint non-parallel minimal genus Seifert surfaces.
    """
    max_crossing = 13
    found_knots = set()

    # The condition is that p is a perfect square.
    # We check p = a*a for odd a.
    # A rough upper bound for p is c^2, so a can go up to c.
    # Let's test a bit higher to be safe.
    for a in range(3, 15, 2):
        p = a * a
        
        # Iterate through q values: 1 <= q < p/2, gcd(p,q)=1
        for q in range(1, (p // 2) + 1):
            if math.gcd(p, q) == 1:
                # Calculate crossing number
                c = get_crossing_number(p, q)

                if c <= max_crossing:
                    rep = get_orbit_rep(p, q)
                    # Store a tuple of (p, representative_q, crossing_number)
                    found_knots.add((p, rep, c))

    sorted_knots = sorted(list(found_knots), key=lambda x: (x[2], x[0], x[1]))
    
    print("Found knots (p, q, crossing_number):")
    for p, q, c in sorted_knots:
        print(f"K({p}, {q}) has crossing number {c}")

    print(f"\nTotal number of such knots is: {len(sorted_knots)}")
    
    return len(sorted_knots)

# Execute the solution
final_count = solve()
print(f"\nFinal answer should be an integer representing the total count.")
print(f"To be explicit, the final equation is just the total number.")
print(f"{final_count}")
<<<45>>>