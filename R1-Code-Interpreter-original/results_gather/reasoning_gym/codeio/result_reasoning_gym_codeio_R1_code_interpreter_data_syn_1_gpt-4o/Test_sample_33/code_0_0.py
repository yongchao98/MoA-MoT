from math import isclose

def find_n_k():
    target_keyspacenomin = 100000
    for n in range(1, 100):  # Arbitrarily chosen upper limit for n
        k = round(target_keyspacenomin ** (1/n))
        if isclose(n ** k, target_keyspacenomin, rel_tol=1e-9) and k < n:
            return {"n": n, "k": k}
    return None

result = find_n_k()
print(result)