import collections

def solve():
    """
    Calculates the number of true boolean expressions of length 5.
    """
    # T[n] stores the number of True expressions of length n.
    # F[n] stores the number of False expressions of length n.
    T = collections.defaultdict(int)
    F = collections.defaultdict(int)

    # Base case: n = 1
    # Expressions: "T", "F"
    T[1] = 1  # "T"
    F[1] = 1  # "F"

    print("Calculating counts for n=1 to 5...")
    print("n=1: T[1]=1, F[1]=1")

    # Dynamic programming for n from 2 to 5
    for n in range(2, 6):
        # 1. Expressions of the form !E, where len(E) = n-1
        t_from_not = F[n-1]
        f_from_not = T[n-1]

        # 2. Expressions of the form (E), where len(E) = n-2
        # This is only possible for n >= 3, since min len(E) is 1.
        t_from_paren = 0
        f_from_paren = 0
        if n >= 3:
            t_from_paren = T[n-2]
            f_from_paren = F[n-2]

        # 3. Expressions of the form E1 op E2, where len(E1)+len(E2) = n-1
        t_from_ops = 0
        f_from_ops = 0
        # Iterate over possible lengths k for the left expression E1
        for k in range(1, n - 1):
            len_e1 = k
            len_e2 = n - 1 - k

            # Total expressions of these lengths
            c1 = T[len_e1] + F[len_e1]
            c2 = T[len_e2] + F[len_e2]
            total_combinations = c1 * c2
            
            # Count for '&' operator
            t_and = T[len_e1] * T[len_e2]
            f_and = total_combinations - t_and
            t_from_ops += t_and
            f_from_ops += f_and

            # Count for '|' operator
            f_or = F[len_e1] * F[len_e2]
            t_or = total_combinations - f_or
            t_from_ops += t_or
            f_from_ops += f_or
            
        T[n] = t_from_not + t_from_paren + t_from_ops
        F[n] = f_from_not + f_from_paren + f_from_ops
        print(f"n={n}: T[{n}]={T[n]}, F[{n}]={F[n]}")
        
        if n == 5:
            print("\n--- Breakdown for T[5] ---")
            print(f"1. From expressions like !E (where E has length 4):")
            print(f"   The number of false expressions of length 4 is F[4] = {F[4]}. So, we have {t_from_not} true expressions.")
            
            print(f"\n2. From expressions like (E) (where E has length 3):")
            print(f"   The number of true expressions of length 3 is T[3] = {T[3]}. So, we have {t_from_paren} true expressions.")
            
            print(f"\n3. From expressions like E1 op E2 (where len(E1)+len(E2)=4):")
            # Decomposing t_from_ops for explanation
            t_k13, t_k22, t_k31 = 0, 0, 0
            # k=1, len_e2=3
            c1, c3 = T[1] + F[1], T[3] + F[3]
            t_k13 += T[1] * T[3] # for &
            t_k13 += c1*c3 - F[1]*F[3] # for |
            # k=2, len_e2=2
            c2 = T[2] + F[2]
            t_k22 += T[2] * T[2] # for &
            t_k22 += c2*c2 - F[2]*F[2] # for |
            # k=3, len_e2=1
            t_k31 = t_k13 # Symmetric
            
            print(f"   - len(E1)=1, len(E2)=3: {t_k13} expressions")
            print(f"   - len(E1)=2, len(E2)=2: {t_k22} expressions")
            print(f"   - len(E1)=3, len(E2)=1: {t_k31} expressions")
            print(f"   Total from operators: {t_k13} + {t_k22} + {t_k31} = {t_from_ops}")

            print("\n--- Final Calculation ---")
            final_sum_str = f"T[5] = (from !) + (from ()) + (from ops)\n"
            final_sum_str += f"T[5] = {t_from_not} + {t_from_paren} + {t_from_ops}\n"
            final_sum_str += f"T[5] = {T[5]}"
            print(final_sum_str)
            
            final_answer = T[5]
            return final_answer

final_answer = solve()
print(f"<<<{final_answer}>>>")