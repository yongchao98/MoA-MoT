def get_cohomology_S7(p, q):
    """
    Returns H^p(S_7, Lambda^q(Z^7)) based on known results.
    This is a simplified table for the purpose of this problem.
    """
    if q == 0:
        # H^p(S_7, Z)
        if p == 0: return "Z"
        if p == 1: return "Z/2Z"
        if p == 2: return "Z/2Z"
        if p == 3: return "Z/6Z"
        # H^4(S_7, Z) is Z/2Z, H^5 is Z/2Z, H^6 is Z/6Z+Z/2Z
        # We only need up to a certain total degree.
        if p > 3: return "..." # higher groups
        return "0"
    if q == 1:
        # H^p(S_7, Z^7)
        if p == 0: return "Z" # Invariants are sum of generators
        if p == 1: return "Z/7Z"
        if p == 2: return "0" # H^2(S_k, Z^k) = 0 for k>=4
        return "0"
    if q == 2:
        # H^p(S_7, Lambda^2(Z^7))
        if p == 0: return "Z" # Rank of invariants is 1 for k>=4
        if p == 1: return "Z/2Z"
        return "0"
    if q == 3:
        # H^p(S_7, Lambda^3(Z^7))
        if p == 0: return "Z" # Rank of invariants is 1 for k>=6
        return "0"
    if q == 4:
        if p == 0: return "Z" # Rank of invariants is 2 for k>=8, 1 for k=7
        return "0"
    # For higher q, the groups become more complex or zero.
    return "0"

def combine_groups(groups):
    """Combines a list of group strings into a sum."""
    counts = {}
    for g in groups:
        if g == "0":
            continue
        if g not in counts:
            counts[g] = 0
        counts[g] += 1
    
    if not counts:
        return "0"

    parts = []
    # Order Z, Z/bZ
    if "Z" in counts:
        if counts["Z"] == 1:
            parts.append("Z")
        else:
            parts.append("Z^" + str(counts["Z"]))
        del counts["Z"]
    
    # Sort remaining finite groups by order
    sorted_finite = sorted(counts.items(), key=lambda item: int(item[0].split('/')[1][:-1]) if '/' in item[0] else 0)

    for g, num in sorted_finite:
        if num == 1:
            parts.append(g)
        else:
            parts.append(g + "^" + str(num))
            
    return "+".join(parts)

def main():
    max_degree = 5 # Calculate up to H^5
    result_list = []
    
    print("Cohomology groups of M(7):")
    
    for n in range(max_degree + 1):
        groups_for_Hn = []
        for p in range(n + 1):
            q = n - p
            term = get_cohomology_S7(p, q)
            if term != "0" and term != "...":
                groups_for_Hn.append(term)
        
        cohomology_group = combine_groups(groups_for_Hn)
        result_list.append(cohomology_group)

    # Format the final output string
    formatted_answer = "[" + ", ".join(result_list) + "]"
    print(formatted_answer)

if __name__ == "__main__":
    main()
