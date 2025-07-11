import math

def calculate_h11(r, l, k_C):
    """Calculates the Hodge number h^{1,1}(M)."""
    return r + 1 + l * k_C

def solve_hodge_number():
    """
    Finds the maximal Hodge number h^{1,1}(M) by searching through known
    configurations of non-symplectic involutions on K3 surfaces.
    """
    # For a curve C of genus 2, the maximal number of fixed points
    # for an involution is 6 (hyperelliptic involution).
    k_C = 6
    print(f"Maximizing h^{1,1}(M) = r + 1 + l * k_C")
    print(f"Choosing the hyperelliptic involution on C gives k_C = {k_C} fixed points.")
    print("-" * 30)

    # Known configurations for the fixed locus S^rho = C_g U k*P^1
    # for a non-symplectic involution on a K3 surface S.
    # The format is r: [(g_1, k_1), (g_2, k_2), ...],
    # where g-k = 11-r must be satisfied.
    configurations = {
        10: [(1, 0)],
        8:  [(3, 0), (4, 1)],
        6:  [(5, 0), (6, 1), (7, 2)],
        4:  [(7, 0), (8, 1), (9, 2), (10, 3)],
        2:  [(9, 0), (10, 1), (11, 2), (12, 3), (13, 4)],
    }

    max_h11 = 0
    best_config = {}

    print("Searching through known configurations (r, g, k):")
    # r is the rank of the invariant part of H^2(S, Z)
    for r in configurations:
        # g is the genus of the main component of the fixed locus S^rho
        # k is the number of rational curves in S^rho
        for g, k in configurations[r]:
            # Verify the configuration with the constraint g - k = 11 - r
            if g - k != 11 - r:
                continue

            # l is the number of connected components of the fixed locus S^rho
            l = k + 1
            h11 = calculate_h11(r, l, k_C)

            if h11 > max_h11:
                max_h11 = h11
                best_config = {'r': r, 'g': g, 'k': k, 'l': l, 'k_C': k_C}

    print("-" * 30)
    print("The optimal configuration found is:")
    print(f"  Rank of invariant lattice r = {best_config['r']}")
    print(f"  Fixed locus S^rho has a curve of genus g = {best_config['g']} and k = {best_config['k']} rational curves.")
    print(f"  This gives l = k + 1 = {best_config['l']} connected components for S^rho.")
    
    r_final = best_config['r']
    l_final = best_config['l']
    k_C_final = best_config['k_C']
    h11_final = best_config['r'] + 1 + best_config['l'] * best_config['k_C']

    print("\nThe maximal value is calculated as follows:")
    print(f"h^{{1,1}}(M) = r + 1 + l * k_C = {r_final} + 1 + {l_final} * {k_C_final} = {h11_final}")
    
    return h11_final

if __name__ == '__main__':
    max_val = solve_hodge_number()
    # The final answer is wrapped according to the required format
    # print(f"\n<<<{max_val}>>>")