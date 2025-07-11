import math

def solve_parliament_k():
    """
    Calculates the maximum integer value of K for the parliament design.
    """
    # 1. Key Parameters
    num_members = 791
    num_sections = 61
    N = math.ceil(num_members / num_sections)  # Number of rows in the worst-case section

    r1 = 3.0  # Radius of the first row
    d = 1.5   # Depth per row
    
    r_N = r1 + (N - 1) * d

    # 2. Find the binding constraint on K
    # The analysis of the worst-case line of sight (observer in row N, opposite section,
    # to speaker in row 1) being obstructed by a person in row j of the speaker's
    # section leads to a boundary condition for K.
    # Let's check the constraint for every possible obstructing row j.
    
    # The condition K <= f(j) must hold for all j.
    # Therefore, K must be <= min(f(j) for j in 2..N)
    
    min_K_bound = float('inf')
    
    print("Deriving the boundary value for K from each potential obstruction...")
    
    # We loop through each potential obstructing person in the speaker's section
    for j in range(2, N + 1):
        r_j = r1 + (j - 1) * d
        
        # The derivation results in the following boundary for K for each row j
        # The derivation relies on the line of sight from the observer at -r_N
        # just clearing the head of a person seated at +r_j.
        # This gives K <= 2 * (r_j + r_N) / ( 1/(r_N+r1) + 1/(r_j-r1) )
        numerator = 2 * (r_j + r_N)
        denominator = (1 / (r_N + r1)) + (1 / (r_j - r1))
        K_bound = numerator / denominator
        
        print(f"For a person in row j={j} (r={r_j}m) to not obstruct, K must be <= {K_bound:.2f}")

        if K_bound < min_K_bound:
            min_K_bound = K_bound
            
    # The final value for K must satisfy the tightest constraint, which is the minimum of all calculated bounds.
    max_K_integer = math.floor(min_K_bound)

    print("\nTo ensure visibility for all seated members, K must be less than or equal to the minimum of these bounds.")
    print(f"The tightest constraint is K <= {min_K_bound:.2f}")
    print(f"The maximum integer value K can take is: {max_K_integer}")
    print("\nFinal Answer Equation:")
    print("N = ceil(791 / 61) = 13")
    print("r_1 = 3m")
    print("d = 1.5m")
    print(f"r_N = r_1 + (N-1)*d = 3 + (13-1)*1.5 = {r_N}")
    j_critical = 2
    r_critical = r1 + (j_critical - 1) * d
    print(f"The critical obstruction is at row j={j_critical} (r={r_critical}m).")
    print(f"max_K = floor(2 * (r_j + r_N) / (1/(r_N+r_1) + 1/(r_j-r_1)))")
    print(f"max_K = floor(2 * ({r_critical} + {r_N}) / (1/({r_N}+{r1}) + 1/({r_critical}-{r1})))")
    print(f"max_K = floor({min_K_bound})")
    
    
    print(f"<<<{max_K_integer}>>>")

solve_parliament_k()