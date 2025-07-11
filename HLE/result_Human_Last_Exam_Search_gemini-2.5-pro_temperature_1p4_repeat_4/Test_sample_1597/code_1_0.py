def get_prime_factor_sets(n):
    factor_sets = [set() for _ in range(n + 1)]
    for i in range(2, n + 1):
        if not factor_sets[i]:
            for j in range(i, n + 1, i):
                factor_sets[j].add(i)
    return factor_sets

limit = 1000
prime_factor_sets = get_prime_factor_sets(limit)

allowed_pairs_count = 0
for a in range(1, limit + 1):
    for b in range(1, limit + 1):
        supp_a = prime_factor_sets[a]
        supp_b = prime_factor_sets[b]
        if supp_a.issubset(supp_b) or supp_b.issubset(supp_a):
            allowed_pairs_count += 1

print(allowed_pairs_count)