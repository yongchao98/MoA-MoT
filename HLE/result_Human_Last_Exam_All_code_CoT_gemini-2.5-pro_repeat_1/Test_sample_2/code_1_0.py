import math

# A simple class to represent finitely generated abelian groups like Z, Z_n, and their direct sums.
class AbelianGroup:
    def __init__(self, free_rank, torsion_orders):
        self.free_rank = free_rank
        # Torsion part stored as a list of orders, e.g., [2, 6] for Z_2 + Z_6
        self.torsion_orders = sorted([o for o in torsion_orders if o > 1])

    def __repr__(self):
        parts = []
        if self.free_rank > 0:
            if self.free_rank == 1:
                parts.append("Z")
            else:
                parts.append(f"Z^{self.free_rank}")
        for order in self.torsion_orders:
            parts.append(f"Z_{order}")
        if not parts:
            return "0"
        return " + ".join(parts)

    def __add__(self, other):
        return AbelianGroup(self.free_rank + other.free_rank, self.torsion_orders + other.torsion_orders)

Z = AbelianGroup(1, [])
Z2 = AbelianGroup(0, [2])
Z6 = AbelianGroup(0, [6])

# Helper function for tensor product of two abelian groups.
def tensor(g1, g2):
    free = g1.free_rank * g2.free_rank
    torsion = []
    # Z_n tensor Z_m = Z_gcd(n,m)
    for t1 in g1.torsion_orders:
        for t2 in g2.torsion_orders:
            torsion.append(math.gcd(t1, t2))
    # Z tensor A = A
    free += g1.free_rank * g2.free_rank
    new_free_rank = g1.free_rank * g2.free_rank
    
    new_torsion_orders = []
    
    # Z x A -> A
    for _ in range(g1.free_rank):
        new_torsion_orders.extend(g2.torsion_orders)
    for _ in range(g2.free_rank):
        new_torsion_orders.extend(g1.torsion_orders)

    # Z_n x Z_m -> Z_gcd(n,m)
    for t1 in g1.torsion_orders:
        for t2 in g2.torsion_orders:
            new_torsion_orders.append(math.gcd(t1, t2))
            
    return AbelianGroup(new_free_rank, new_torsion_orders)

# Helper function for Tor functor. Tor(Z_n, Z_m) = Z_gcd(n,m). Tor(Z, A) = 0.
def tor(g1, g2):
    torsion = []
    for t1 in g1.torsion_orders:
        for t2 in g2.torsion_orders:
            torsion.append(math.gcd(t1, t2))
    return AbelianGroup(0, torsion)

# Universal Coefficient Theorem: H_p(X; A) = H_p(X) tensor A + Tor(H_{p-1}(X), A)
def homology_with_coeffs(p, Hp, Hpm1, A):
    term1 = tensor(Hp.get(p, AbelianGroup(0, [])), A)
    term2 = tor(Hp.get(p - 1, AbelianGroup(0, [])), A)
    return term1 + term2

def compute_spin_bordism():
    """
    Computes the reduced 12-th dimensional Spin bordism of BG2.
    """
    # Step 1 & 2: Define input data
    # Spin bordism groups of a point Omega_q^Spin for q <= 12
    Omega_spin = {
        0: Z,
        1: Z2,
        2: Z2,
        3: AbelianGroup(0, []),
        4: Z,
        5: AbelianGroup(0, []),
        6: AbelianGroup(0, []),
        7: AbelianGroup(0, []),
        8: Z + Z,
        9: Z2 + Z2,
        10: Z2,
        11: AbelianGroup(0, []),
        12: Z + Z,
    }

    # Integral homology groups H_p(BG2; Z) for p <= 12
    # Derived from H*(BG2; Z) = Z[y4, y8, y12]/(2y8-y4^2, y4y8-3y12, y8^2-4y4y12)
    # H^4=Z, H^8=Z, H^12=Z+Z_6. UCT gives H_p.
    H_p_BG2 = {
        0: Z,
        4: Z,
        8: Z,
        11: Z6,
        12: Z,
    }
    
    print("Computing the reduced 12-th Spin bordism of BG2, denoted Omega_12^Spin(BG2).")
    print("Using the Atiyah-Hirzebruch Spectral Sequence (AHSS).")
    print("E^2_{p,q} = H_p(BG2; Omega_q^Spin), converging to Omega_{p+q}^Spin(BG2).\n")

    total_group = AbelianGroup(0, [])
    reduced_group = AbelianGroup(0, [])
    
    # Step 3: Calculate E^2_{p,q} for p+q=12
    print("Calculating E^2_{p,q} terms for p+q = 12:")
    for p in range(13):
        q = 12 - p
        if q < 0: continue
        
        Omega_q = Omega_spin.get(q, AbelianGroup(0,[]))
        
        # E^2_{p,q} = H_p(BG2; Omega_q^Spin)
        E2_pq = homology_with_coeffs(p, H_p_BG2, H_p_BG2, Omega_q)

        if E2_pq.free_rank > 0 or E2_pq.torsion_orders:
            print(f"E^2_{p},{q} = H_{p}(BG2; Omega_{q}^Spin) = H_{p}(BG2; {Omega_q}) = {E2_pq}")
            total_group += E2_pq
            if p > 0:
                reduced_group += E2_pq

    # Step 4, 5, 6, 7
    print("\nThe AHSS is known to collapse for BG2, so E^2 = E^infinity.")
    print("The group Omega_12^Spin(BG2) is the direct sum of these E^2 terms.")
    print(f"So, Omega_12^Spin(BG2) = {total_group}")
    
    print("\nThe 'reduced' group is the full group quotiented by the p=0 term.")
    print(f"The term for p=0 is E^2_0,12 = {homology_with_coeffs(0, H_p_BG2, H_p_BG2, Omega_spin[12])}")
    print(f"So, the reduced group is {total_group} / {homology_with_coeffs(0, H_p_BG2, H_p_BG2, Omega_spin[12])}")
    print("\nFinal Result:")
    print(f"The reduced 12-th Spin bordism group of BG2 is: {reduced_group}")
    
    # Final answer in the required format
    # The final answer is Z^4 + Z_2. This will be wrapped by <<<>>>.
    final_answer_str = repr(reduced_group)
    return final_answer_str

final_result = compute_spin_bordism()
# The required output format is just the answer content.
# The code produces "Z^4 + Z_2"
# I will output it as required
final_answer_for_user = "<<<" + final_result + ">>>"
#print(final_answer_for_user)
# The user wants to see the python script, not the direct output from my execution of it.
# So I should not print final_answer_for_user. Just the code block.