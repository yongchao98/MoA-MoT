import math
from collections import defaultdict

def get_pq(cf):
    """Computes the fraction p/q from a continued fraction list."""
    p, q = 1, 0
    for c in reversed(cf):
        p, q = c * p + q, p
    return p, q

def compositions(n, k):
    """Generates all compositions of integer n into k parts."""
    if k == 1:
        yield [n]
        return
    for i in range(1, n - k + 2):
        for rest in compositions(n - i, k - 1):
            yield [i] + rest

def get_mod_inverse(a, m):
    """Computes modular inverse of a modulo m."""
    if math.gcd(a, m) != 1:
        return None
    return pow(a, -1, m)

def count_knots():
    """
    Counts the number of 2-bridge knots with c<=13 having two disjoint
    minimal genus Seifert surfaces, treating mirror images as indistinct.
    """
    max_crossing = 13
    knot_fractions = defaultdict(set)

    for c in range(4, max_crossing + 1, 2):
        s = c // 2
        for k in range(2, c + 1, 2):
            if k > s:
                continue
            for comp in compositions(s, k):
                cf = [2 * a for a in comp]
                p, q = get_pq(cf)
                if p % 2 != 0 and math.gcd(p, q) == 1:
                    knot_fractions[p].add(q)
    
    total_knots = 0
    knot_counts = {}

    # Sort p for deterministic output
    sorted_p = sorted(knot_fractions.keys())

    all_p_q_classes = defaultdict(list)
    
    for p in sorted_p:
        qs = sorted(list(knot_fractions[p]))
        
        # Determine crossing number for this p
        # Find one q and calculate its CF to get crossing number
        # This is inefficient, but clear. For a given p, all q will have same c
        q_sample = qs[0]
        temp_p, temp_q = p, q_sample
        cf = []
        while temp_q != 0:
            c_i = temp_p // temp_q
            cf.append(c_i)
            temp_p, temp_q = temp_q, temp_p % temp_q
        c = sum(cf)
        
        # Partition qs into equivalence classes
        classes = []
        while qs:
            q_rep = qs.pop(0)
            
            q_inv = get_mod_inverse(q_rep, p)
            if q_inv is None: # Should not happen as gcd(p,q)=1
                continue

            orbit = {q_rep, q_inv, (p - q_rep) % p, (p - q_inv) % p}
            
            new_class = {q_rep}
            
            # Find all other q's in the same orbit
            remaining_qs = []
            for q_item in qs:
                if q_item in orbit:
                    new_class.add(q_item)
                else:
                    remaining_qs.append(q_item)
            qs = remaining_qs
            classes.append(new_class)
        
        all_p_q_classes[c].append(len(classes))
    
    print("Number of such knots found for each crossing number:")
    final_sum_str = []
    
    # Sort by crossing number
    for c in sorted(all_p_q_classes.keys()):
        count = sum(all_p_q_classes[c])
        knot_counts[c] = count
        total_knots += count
        print(f"Crossing number {c}: {count}")
        final_sum_str.append(str(count))

    print("\nTotal number of knots:")
    print(f"= {' + '.join(final_sum_str)}")
    print(f"= {total_knots}")

if __name__ == '__main__':
    count_knots()
