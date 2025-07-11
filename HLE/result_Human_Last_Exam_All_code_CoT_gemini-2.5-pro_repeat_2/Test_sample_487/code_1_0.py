import math

def calculate_h11_M():
    """
    Calculates the maximal value of the Hodge number h^{1,1}(M).
    
    M is a smooth resolution of the quotient (S x C) / (rho x psi), where:
    S is a K3 surface.
    C is a complex curve of genus 2.
    rho is a non-symplectic involution on S.
    psi is an involution on C.
    """

    # Step 1: Formulate h^{1,1}(M)
    # h^{1,1}(M) = h^{1,1}(quotient) + N_res
    # where N_res is the number of exceptional divisors from the resolution.

    # Step 2: Calculate h^{1,1}(quotient)
    # By Künneth formula, H^{1,1}(S x C) = (H^{1,1}(S) tensor H^{0,0}(C)) + (H^{0,0}(S) tensor H^{1,1}(C)).
    # h^{1,1}(S) = 20, h^{1,1}(C) = 1.
    # The invariant part under G = rho x psi is h^{1,1}_+(rho) + 1.
    # h^{1,1}(quotient) = h^{1,1}_+(rho) + 1

    # Step 3: Calculate N_res
    # The singular locus is the image of Fix(rho) x Fix(psi).
    # N_res = (number of components of Fix(rho)) * (number of fixed points of psi)
    # Let k = number of components of Fix(rho)
    # Let N_psi = number of fixed points of psi
    # N_res = k * N_psi

    # Step 4: Relate h^{1,1}_+(rho) to the geometry of Fix(rho)
    # The Lefschetz fixed-point formula for rho on S gives:
    # chi(Fix(rho)) = 2 * h^{1,1}_+(rho) - 20
    # Also, chi(Fix(rho)) = sum(2 - 2*g_i) for each component C_i of Fix(rho) with genus g_i.
    # Let sum_g = sum of genera of components of Fix(rho).
    # 2*k - 2*sum_g = 2 * h^{1,1}_+(rho) - 20
    # h^{1,1}_+(rho) = 10 + k - sum_g

    # Step 5: Combine everything to get the formula for h^{1,1}(M)
    # h^{1,1}(M) = (10 + k - sum_g) + 1 + (k * N_psi)
    # h^{1,1}(M) = 11 + k * (1 + N_psi) - sum_g
    print("The formula for h^{1,1}(M) is derived as: h^{1,1}(M) = 11 + k * (1 + N_psi) - sum_g")
    print("where k is the number of components of Fix(rho), sum_g is the sum of their genera, and N_psi is the number of fixed points of psi.\n")

    # Step 6: Maximize the expression by choosing psi and rho.

    # Choice for psi on a genus 2 curve C:
    # An involution can have 2 or 6 fixed points (N_psi).
    # To maximize h^{1,1}(M), we need to maximize N_psi.
    N_psi = 6  # Hyperelliptic involution
    print(f"To maximize, we choose the involution psi on the curve C to be the hyperelliptic one, which gives the maximal number of fixed points.")
    print(f"Maximal N_psi = {N_psi}\n")
    
    # The formula becomes: h^{1,1}(M) = 11 + k * (1 + 6) - sum_g = 11 + 7*k - sum_g
    # Let's call the term to maximize T(rho) = 7*k - sum_g

    # Choice for rho on a K3 surface S:
    # We analyze the cases for Fix(rho) from Nikulin's classification of non-symplectic involutions.
    print("To maximize further, we analyze the possible fixed loci Fix(rho) for a non-symplectic involution on a K3 surface:")
    
    cases = []
    
    # Case 1: Fix(rho) is empty.
    k1 = 0
    sum_g1 = 0
    term1 = 7 * k1 - sum_g1
    h11_1 = 11 + term1
    cases.append({'name': 'Empty set', 'k': k1, 'sum_g': sum_g1, 'term': term1, 'h11': h11_1})

    # Case 2: Fix(rho) is a single smooth curve of genus g, where 0 <= g <= 10.
    # To maximize 7*k - sum_g = 7*1 - g = 7 - g, we must minimize g.
    k2 = 1
    g2 = 0
    sum_g2 = g2
    term2 = 7 * k2 - sum_g2
    h11_2 = 11 + term2
    cases.append({'name': 'Single rational curve (g=0)', 'k': k2, 'sum_g': sum_g2, 'term': term2, 'h11': h11_2})
    
    # Case 3: Fix(rho) is the disjoint union of two smooth rational curves.
    k3 = 2
    sum_g3 = 0 + 0
    term3 = 7 * k3 - sum_g3
    h11_3 = 11 + term3
    cases.append({'name': 'Two rational curves', 'k': k3, 'sum_g': sum_g3, 'term': term3, 'h11': h11_3})

    max_h11 = 0
    best_case = None
    
    for case in cases:
        print(f"- Case: Fix(rho) is a {case['name']}.")
        print(f"  k = {case['k']}, sum_g = {case['sum_g']}")
        print(f"  Maximizing term = 7 * {case['k']} - {case['sum_g']} = {case['term']}")
        print(f"  Resulting h^{1,1}(M) = 11 + {case['term']} = {case['h11']}\n")
        if case['h11'] > max_h11:
            max_h11 = case['h11']
            best_case = case
            
    print("Comparing the results, the maximum value is obtained in the last case.")
    final_term = best_case['term']
    print(f"The maximal value of the term (7k - sum_g) is {final_term}.")
    print(f"The maximal value of h^(1,1)(M) is therefore 11 + {final_term} = {max_h11}.")

    return max_h11

# Run the calculation and print the final answer.
max_value = calculate_h11_M()
print(f"\nFinal calculation is: 11 + (7 * 2 - 0) = {max_value}")
print(f"<<<{max_value}>>>")

calculate_h11_M()