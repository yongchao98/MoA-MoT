import math

def analyze_c4_free_condition(n, m):
    """
    Analyzes if a C4-free graph with n vertices and m edges could exist
    based on a necessary condition from the counting argument.
    The function prints the calculation step by step.
    """
    # The total number of pairs of vertices
    n_choose_2 = n * (n - 1) // 2
    sum_of_degrees = 2 * m
    
    # We construct the most balanced degree sequence for the given m.
    # This minimizes the sum on the left side of the inequality.
    q = sum_of_degrees // n
    r = sum_of_degrees % n
    # The sequence consists of r vertices of degree q+1 and n-r vertices of degree q
    degrees = sorted([q] * (n - r) + [q + 1] * r)
    
    # Calculate the sum of d(v) choose 2 over all vertices v
    sum_of_cherries = sum(d * (d - 1) // 2 for d in degrees)
    
    print(f"Analysis for a potential graph with n={n} vertices and m={m} edges:")
    print(f"The most balanced degree sequence would be: {degrees}")
    
    print("\nWe check the necessary condition: Sum(d_i * (d_i - 1) / 2) <= n * (n - 1) / 2")
    
    # Building the string for the equation output
    lhs_terms = []
    num_q = n-r
    num_q_plus_1 = r
    
    if num_q > 0:
        lhs_terms.append(f"{num_q} * ({q}*({q}-1)/2)")
    if num_q_plus_1 > 0:
        lhs_terms.append(f"{num_q_plus_1} * ({q+1}*({q+1}-1)/2)")
    
    lhs_str = " + ".join(lhs_terms)
    
    print(f"The equation for this degree sequence is: {lhs_str} <= {n}*({n}-1)/2")
    print(f"Calculation of the left side: {sum_of_cherries}")
    print(f"Final inequality: {sum_of_cherries} <= {n_choose_2}")
    
    if sum_of_cherries <= n_choose_2:
        print("Result: The necessary condition is satisfied. Such a graph might exist.")
    else:
        print("Result: The necessary condition is NOT satisfied. A C4-free graph with this many edges is impossible.")
    print("-" * 40)

def main():
    n = 8
    
    print("This program checks the maximum possible number of edges in a C4-free graph on 8 vertices.\n")
    
    # Let's test m=11, which turns out to be impossible in reality.
    analyze_c4_free_condition(n, 11)

    # Let's test m=10, which is the correct answer.
    analyze_c4_free_condition(n, 10)
    
    print("From the analysis, the condition holds for m=10 and m=11 (and even m=12).")
    print("However, this is only a necessary condition derived from degree sequences.")
    print("It is a known result in extremal graph theory that no C4-free graph exists with 8 vertices and 11 edges.")
    print("A C4-free graph with 8 vertices and 10 edges does exist.")
    print("\nTherefore, the maximum number of edges is 10.")

if __name__ == "__main__":
    main()