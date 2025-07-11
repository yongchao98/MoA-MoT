import math

def solve_question_1():
    """
    Solves: Is it true that for a d-neighbor N of Z^12, R2(M) can be of type A_11?
    """
    n = 12
    # R2(Z^n) is of type D_n. We want R2(M) to be A_{11}.
    # The roots of A_{11} in Z^{12} can be taken as {e_i - e_j | i != j}.
    # The roots of D_{12} are {+/- e_i +/- e_j | i < j}.
    # We want a sublattice M that contains {e_i - e_j} but not {e_i + e_j}.
    # Let M = {x in Z^n | x . u == 0 (mod d)}.
    # 1. (e_i - e_j) . u == 0 (mod d)  => u_i - u_j == 0 => u_i == u_j (mod d).
    #    Let's say u_i == c (mod d) for all i.
    # 2. (e_i + e_j) . u != 0 (mod d) => u_i + u_j != 0 => 2c != 0 (mod d).
    # 3. For a d-neighbor to exist, we need d | u . u for some primitive lift of u.
    #    Let's test u = (1, 1, ..., 1). u is primitive. u_i = 1 for all i, so we can set c=1.
    #    u . u = n = 12. We need d | 12.
    #    Condition 2 becomes 2*1 != 0 (mod d), so d cannot be 1 or 2.
    #    Let's test d=3. d|12 -> yes. d != 1,2 -> yes.
    #    So d=3 with u=(1,...,1) is a valid choice.
    #    This shows a neighbor can be constructed where R2(M) is of type A_11.
    answer = "Yes"
    print("Part (a):")
    print("For n=12, we want to know if the visible root system R2(M) can be A_11.")
    print("The root system of Z^12 is D_12, which consists of vectors (+/- e_i +/- e_j).")
    print("The root system A_11 can be represented by vectors {e_i - e_j}.")
    print("A_11 is a subsystem of D_12. We need to check if there is a sublattice M = Z^12 intersect N such that R2(M) is precisely A_11.")
    print("Let's construct M as M = {x in Z^12 | x . u == 0 (mod d)} for a vector u and index d.")
    print("For all roots of A_11 to be in M, we need (e_i - e_j) . u == 0 (mod d), which means u_i == u_j (mod d).")
    print("To exclude other roots like e_i + e_j, we need (e_i + e_j) . u != 0 (mod d), which means u_i + u_j != 0 (mod d).")
    print("Let's try the simple primitive vector u = (1, 1, ..., 1). Then u_i = 1 for all i.")
    print("The first condition holds. The second becomes 1 + 1 = 2 != 0 (mod d), so d cannot be 2.")
    print("A d-neighbor can be constructed if d divides u.u. Here u.u = 12.")
    print(f"We need to find a d such that d > 2 and d divides 12. Let's choose d = 3.")
    print(f"With d=3 and u=(1,...,1), M contains all e_i-e_j (since (e_i-e_j).u = 1-1=0) but not e_i+e_j (since (e_i+e_j).u = 1+1=2 != 0 mod 3).")
    print("Therefore, it is possible.")
    return answer

def solve_question_2():
    """
    Solves: Can the visible root system R2(M) of a d-neighbor N of Z^15 contain a D_7 component?
    """
    n = 15
    k = 7 # for D_7
    # We want R2(M) to contain a D_k component, where k=7.
    # This means all roots {+/- e_i +/- e_j | 1 <= i < j <= k} must be in M.
    # Let M = {x in Z^n | x . u == 0 (mod d)}.
    # 1. (+/- e_i +/- e_j) . u == 0 (mod d) for i,j in {1..k}.
    #    This implies u_i == u_j (mod d) and 2*u_i == 0 (mod d) for i,j in {1..k}.
    # To have a D_k *component*, it must be separable from other coordinates {k+1...n}.
    # Let's try d=2.
    # Then 2*u_i == 0 (mod 2) is always true.
    # The condition u_i == u_j (mod 2) means u_1, ..., u_k must have the same parity.
    # Separation from j>k requires u_i and u_j must have different parities.
    # We also need d | u.u, so u.u must be even.
    # k * parity_1 + (n-k) * parity_2 == 0 (mod 2).
    # n=15, k=7. So 7 * p1 + 8 * p2 == 0 (mod 2).
    # This forces p1=0.
    # So u_1..u_7 must be even, and u_8..u_15 must be odd.
    # Let u = (2, 0..0, 1, 1..1). k=7, n=15. u has u_1=2, u_2..u_7=0, u_8..u_15=1.
    # u is primitive.
    # u.u = 2^2 + 6*0^2 + 8*1^2 = 4 + 8 = 12.
    # We need d | u.u. With d=2, 2|12. Yes.
    answer = "yes"
    print("\nPart (b):")
    print(f"For n=15, we want to know if R2(M) can contain a D_7 component.")
    print("For a D_k component on coordinates {1..k}, the construction vector u must satisfy:")
    print("1. u_i == u_j (mod d) and 2*u_i == 0 (mod d) for all i,j in {1..k}.")
    print("2. To be a separate component, roots mixing coordinates inside and outside {1..k} must be excluded.")
    print("3. For the neighbor to exist, d must divide u.u for some primitive integer lift of u.")
    print("Let's try d=2. Condition (1) becomes: u_1,..,u_7 must have the same parity.")
    print("Separation implies u_1,..,u_7 and u_8,..,u_15 must have different parities.")
    print("The neighbor condition (d | u.u) means u.u must be even. This forces the parity of u_1..u_7 to be even.")
    print("So we need u_1..u_7 to be even and u_8..u_15 to be odd.")
    print(f"Let's choose a primitive vector u satisfying this: u = (2,0,0,0,0,0,0, 1,1,1,1,1,1,1,1).")
    print(f"The vector u is primitive. u.u = 2^2 + 8*1^2 = 12.")
    print(f"d=2 divides u.u=12. So a 2-neighbor exists.")
    print("With this u, the conditions for a D_7 component on the first 7 coordinates are met.")
    print("Therefore, it is possible.")
    return answer


def solve_question_3():
    """
    Solves: For n = 18 and d = 5, is it possible for R2(M) to include more than one D_n component?
    """
    n = 18
    d = 5
    # To include two D components, say D_k1 on index set I1 and D_k2 on I2 (disjoint).
    # Let M = {x in Z^n | x . u == 0 (mod d)}.
    # For D_k1: u_i == c1 (mod d) and 2*c1 == 0 (mod d) for i in I1.
    # For D_k2: u_j == c2 (mod d) and 2*c2 == 0 (mod d) for j in I2.
    # For them to be separate components, mixed roots must be excluded.
    # This requires c1 +/- c2 to be non-zero (mod d).
    # Let's solve 2c == 0 (mod 5). Since 5 is prime and does not divide 2,
    # we must have c == 0 (mod 5).
    # This means c1=0 and c2=0.
    # Then c1 +/- c2 = 0. This violates the separation condition.
    answer = "no"
    print("\nPart (c):")
    print(f"For n=18, d=5, we want to know if R2(M) can have more than one D_k component.")
    print("Suppose R2(M) contains D_{k1} on coordinate set I_1 and D_{k2} on I_2.")
    print("This requires the construction vector u to satisfy certain congruences modulo d=5.")
    print(f"For a D_k component on a set I, we need u_i == c (mod 5) for all i in I, with 2*c == 0 (mod 5).")
    solutions = [c for c in range(d) if (2 * c) % d == 0]
    print(f"Solving 2*c == 0 (mod 5) for c in {{0,1,2,3,4}} gives c = {solutions[0]}.")
    c = solutions[0]
    print(f"So for the D_{k1} component, we need u_i == {c} (mod 5) for i in I_1.")
    print(f"For the D_{k2} component, we need u_j == {c} (mod 5) for j in I_2.")
    print("For these components to be separate, roots that mix coordinates from I_1 and I_2 must be excluded from M.")
    print(f"Let's check a mixed root, v = e_i + e_j, where i is in I_1 and j is in I_2.")
    print(f"v.u = u_i + u_j (mod 5). Based on our findings, this is {c} + {c} = {c+c} (mod 5).")
    print("Since v.u is 0 mod 5, this mixed root is in M. This connects the two supposed components.")
    print("The root system on I_1 union I_2 is therefore D_{k1+k2}, a single component.")
    print("Therefore, it is not possible for R2(M) to have more than one D-type component.")
    return answer


def main():
    """
    Main function to answer the user's three-part question.
    """
    ans1 = solve_question_1()
    ans2 = solve_question_2()
    ans3 = solve_question_3()
    final_answer = f"(a) [{ans1}]; (b) [{ans2}]; (c) [{ans3}]."
    print("\n" + "="*20)
    print("Final Answer:")
    print(final_answer)


main()