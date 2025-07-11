def is_shifted(family, n_max):
    """
    Checks if a family of sets is shifted.
    A family F is shifted if for any F in F, and any i in F, j not in F with j < i,
    the set (F - {i}) U {j} is also in F.
    """
    family_of_frozensets = {frozenset(s) for s in family}
    for f_set in family:
        for i in f_set:
            for j in range(1, i):
                if j not in f_set:
                    shifted_set = (f_set - {i}) | {j}
                    if frozenset(shifted_set) not in family_of_frozensets:
                        return False
    return True

def is_t_intersecting(family, t):
    """
    Checks if a family is t-intersecting.
    This means for any two sets F1, F2 in the family, |F1 intersect F2| >= t.
    """
    for i in range(len(family)):
        for j in range(i, len(family)):
            f1 = family[i]
            f2 = family[j]
            if len(f1.intersection(f2)) < t:
                return False
    return True

def main():
    """
    Demonstrates the counterexample for question (b).
    """
    # Parameters for the counterexample
    t = 1
    k = 3
    n = 7
    
    print(f"Checking case: t={t}, k={k}, n={n}")
    
    # Check if n >= k + t + 3
    condition_n = k + t + 3
    print(f"Condition is n >= k + t + 3, which is {n} >= {condition_n}.")
    if n < condition_n:
        print("Parameters do not meet the condition.")
        return
    print("Parameters meet the condition.\n")

    # Define the family F
    F_family = [{1, 2, 3}]
    print(f"Counterexample family F = {F_family}")

    # Verify properties of F
    shifted_check = is_shifted(F_family, n)
    intersecting_check = is_t_intersecting(F_family, t + 1)
    
    print(f"Is F shifted? {shifted_check}")
    print(f"Is F {t+1}-intersecting? {intersecting_check}")
    
    if not (shifted_check and intersecting_check):
        print("The family does not satisfy the premises.")
        return
    
    print("\nThe family is a valid shifted {t+1}-intersecting family.\n")
    
    # Compute F^(n)
    F_n = [f_set for f_set in F_family if n not in f_set]
    
    print(f"F^({n}) = {{ F in F : {n} not in F }}")
    print(f"F^({n}) is: {F_n}")
    
    size_F_n = len(F_n)
    print(f"The size |F^({n})| is {size_F_n}.")
    
    # Check the conclusion
    must_be_ge = 3
    print(f"The question is if |F^({n})| must be >= {must_be_ge}.")
    print(f"Is {size_F_n} >= {must_be_ge}? {size_F_n >= must_be_ge}")
    
    if size_F_n < must_be_ge:
        print("\nConclusion: We have found a counterexample. The answer to (b) is No.")

if __name__ == "__main__":
    main()
