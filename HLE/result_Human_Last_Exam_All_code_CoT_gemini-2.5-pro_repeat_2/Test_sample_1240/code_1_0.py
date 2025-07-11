import math

def get_d_component_indices(d):
    """Finds indices j for which S_j can form a D-type component."""
    return [j for j in range(d) if (2 * j) % d == 0]

def solve_q1():
    """Q1: Is it true that for a d-neighbor N of Z^12, R2(M) can be of type A_11?"""
    # We seek a root system of type A_11 for n=12.
    # This can occur if the system is A_{n_k-1}, meaning n_k=12 and all other n_j=0.
    # This requires choosing a characteristic vector w where all components are
    # congruent to k mod d, e.g., w = (k, k, ..., k).
    # For this to be a valid choice, we need gcd(k, d) = 1.
    # To ensure the root system is A_11 and not D_12, we need 2k != 0 (mod d).
    
    # Search for a valid (d, k) pair.
    for d in range(2, 13):
        for k in range(1, d):
            is_primitive = (math.gcd(k, d) == 1)
            is_not_d_type = ((2 * k) % d != 0)
            
            if is_primitive and is_not_d_type:
                # Example: d=3, k=1. w=(1,..,1). gcd(1,3)=1. 2*1 % 3 != 0.
                # This gives partition n_1=12. System is A_{12-1} = A_11.
                return "Yes"
    return "No"

def solve_q2():
    """Q2: Can the visible root system R2(M) of a d-neighbor N of Z^15 contain a D_7 component?"""
    # We need a D_7 component for n=15.
    # This requires a partition S_j with size n_j = 7, where 2j % d == 0.
    # The remaining n - n_j = 15 - 7 = 8 indices must be partitioned into other S_k's.
    # A characteristic vector w for such a partition must exist.
    # Let's assume the 8 remaining indices are in a single set S_k.
    # The existence condition for w is gcd(j, k, d) = 1.
    
    for d in range(2, 16):
        d_indices = get_d_component_indices(d)
        for j in d_indices:
            # We have a candidate j for the D_7 component (n_j = 7).
            # Now, place the remaining 8 elements in S_k.
            for k in range(d):
                if k == j:
                    continue
                
                # Check if this choice of w is valid.
                # gcd(w_1, ..., w_n, d) is gcd(j, k, d).
                if math.gcd(math.gcd(j, k), d) == 1:
                    # Example: d=4, j=2, k=1. 2*2 % 4 == 0. gcd(2, 1, 4) = 1.
                    # R2(M) would be D_7 + A_7. This contains a D_7 component.
                    return "Yes"
    return "No"

def solve_q3():
    """Q3: For n=18 and d=5, is it possible for R2(M) to include more than one D_n component?"""
    d = 5
    # D-type components arise from partitions S_j where 2j is a multiple of d.
    d_indices = get_d_component_indices(d)
    
    # The number of possible D-components is the number of such j's.
    # For d=5, 2j = 0 (mod 5) implies j=0, since gcd(2,5)=1.
    # So there is only one such j.
    if len(d_indices) > 1:
        return "Yes"
    else:
        return "No"

def main():
    """
    Solves the three questions and prints the answer in the specified format.
    """
    q1_ans = solve_q1()
    q2_ans = solve_q2().lower()
    q3_ans = solve_q3().lower()

    final_answer_string = f"(a) [{q1_ans}]; (b) [{q2_ans}]; (c) [{q3_ans}]."
    print(final_answer_string)

if __name__ == "__main__":
    main()