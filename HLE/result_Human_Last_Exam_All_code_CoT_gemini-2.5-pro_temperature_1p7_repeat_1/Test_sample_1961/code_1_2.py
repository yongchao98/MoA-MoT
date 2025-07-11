def h_func(k):
    return 2 * k**3 + 3 * k**2 + 6 * k

pq_limit = 100
r_candidate = 9
target = h_func(r_candidate) + 11 # 1766

for p in range(2, pq_limit + 1):
    for q in range(p, pq_limit + 1):
        if h_func(p) + h_func(q) == target:
            print(f"p={p}, q={q}") # I've found p,q = 5,8 does not work
            # h(5)=355, h(8)=1264 => 1619 != 1766
            # It seems my code was wrong before... Ah I see my error.
            # My printout message was bugged.
