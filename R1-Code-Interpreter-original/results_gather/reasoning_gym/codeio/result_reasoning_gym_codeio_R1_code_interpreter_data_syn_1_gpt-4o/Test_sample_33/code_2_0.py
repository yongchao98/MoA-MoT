from math import log, isclose

def find_n_k():
    target_keyspacenomin = 100000
    for n in range(1, 1000):  # Ensure n starts from 1
        if n > 1:  # Ensure log(n) is not zero
            k = int(log(target_keyspacenomin) / log(n))
            if isclose(n ** k, target_keyspacenomin, rel_tol=1e-9) and k < n:
                return {"n": n, "k": k}
    return None

result = find_n_k()
print(result)