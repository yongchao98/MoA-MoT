import math
from fractions import Fraction

def calculate_c_n(n):
    """
    Calculates the prefactor c_n for a given n.

    The function uses pre-computed data for the number of labeled, biconnected
    graphs on n vertices, partitioned by the number of edges.
    Data is available for n from 2 to 6.
    """
    if not isinstance(n, int) or n < 2:
        print("Error: n must be an integer greater than or equal to 2.")
        return

    # Data source: OEIS A006125
    # T[n][k] is the number of labeled biconnected graphs on n vertices with k edges.
    T = {
        2: {1: 1},
        3: {3: 1},
        4: {4: 3, 5: 6, 6: 1},
        5: {5: 12, 6: 60, 7: 90, 8: 60, 9: 15, 10: 1},
        6: {6: 60, 7: 540, 8: 2220, 9: 5400, 10: 8580, 11: 9108, 12: 6435, 13: 2970, 14: 840, 15: 112},
    }

    if n not in T:
        print(f"Sorry, data for n={n} is not available in this script.")
        return

    print(f"Calculating the prefactor c_{n}:")
    
    # m_max is the number of edges in a complete graph K_n
    m_max = n * (n - 1) // 2
    print(f"For n = {n}, the number of edges in a complete graph K_n is m_max = {m_max}.")

    graph_counts = T[n]
    print("\nThe biconnected graphs on {n} vertices are distributed as follows:")
    for edges, count in graph_counts.items():
        print(f"- {count} graph(s) with {edges} edges.")

    # Calculate the sum C(n)
    C_n_sum = 0
    sum_str_parts = []
    for m, count in graph_counts.items():
        term = count * ((-1)**(m_max - m))
        C_n_sum += term
        sum_str_parts.append(f"{count}*(-1)^({m_max}-{m})")

    sum_str = " + ".join(sum_str_parts)
    print(f"\nThe alternating sum C({n}) = {sum_str}")
    
    # Calculate each term to show the sum
    sum_calc_parts = []
    for m, count in graph_counts.items():
        term_val = count * ((-1)**(m_max - m))
        sum_calc_parts.append(str(term_val))
    sum_calc_str = " + ".join(sum_calc_parts).replace('+ -', '- ')
    
    print(f"C({n}) = {sum_calc_str} = {C_n_sum}")

    # Calculate c_n
    numerator = -(n - 1) * C_n_sum
    denominator = math.factorial(n)
    c_n_frac = Fraction(numerator, denominator)
    
    print(f"\nThe prefactor c_{n} = -({n}-1)/{n}! * C({n})")
    print(f"c_{n} = -{n-1}/{math.factorial(n)} * ({C_n_sum})")
    print(f"c_{n} = {numerator}/{denominator}")
    print(f"c_{n} = {c_n_frac}")

# Example calculation for n=5
calculate_c_n(5)