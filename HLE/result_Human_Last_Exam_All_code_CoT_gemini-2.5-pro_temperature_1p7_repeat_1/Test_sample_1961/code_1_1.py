def h_func(k):
    return 2 * k**3 + 3 * k**2 + 6 * k
r_limit = 100
pq_limit = 100
h_values = {h_func(r): r for r in range(2, r_limit + 1)}
for r_candidate in range(2, r_limit + 1):
    target = h_func(r_candidate) + 11
    found = False
    for p in range(2, pq_limit + 1):
        if found: break
        for q in range(p, pq_limit + 1):
            if h_func(p) + h_func(q) == target:
                # print(f"p={p}, q={q}, r={r_candidate}")
                # print(f"{h_func(p)} + {h_func(q)} = {h_func(p)+h_func(q)}")
                # print(f"{h_func(r_candidate)} + 11 = {h_func(r_candidate)+11}")
                print(r_candidate) # a simpler output
                found = True
                break
    if found:
        break