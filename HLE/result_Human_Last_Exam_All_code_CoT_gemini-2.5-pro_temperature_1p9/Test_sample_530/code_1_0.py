import math

class EndomorphismRing:
    """A simple class to represent an endomorphism ring by its rank."""
    def __init__(self, rank):
        # The rank of the endomorphism ring as a Z-module.
        if not isinstance(rank, int) or rank < 0:
            raise ValueError("Rank must be a non-negative integer.")
        self.rank = rank

class Torus:
    """Represents an algebraic torus T = (G_m)^r."""
    def __init__(self, dimension):
        self.dimension = dimension
        # End(T) is isomorphic to M_r(Z), the ring of r x r integer matrices.
        # The rank of M_r(Z) as a Z-module is r^2.
        self.endomorphism_ring = EndomorphismRing(dimension * dimension)

class AbelianVariety:
    """Represents an abelian variety A."""
    def __init__(self, name, endomorphism_rank):
        self.name = name
        self.endomorphism_ring = EndomorphismRing(endomorphism_rank)

class SemiAbelianVariety:
    """Represents a semi-abelian variety G, as an extension of A by T."""
    def __init__(self, A, T, extension_type, rigid_rank=1):
        self.A = A # underlying abelian variety
        self.T = T # torus
        self.extension_type = extension_type
        self.endomorphism_ring, self.equation = self._calculate_endomorphism_rank(rigid_rank)

    def _calculate_endomorphism_rank(self, rigid_rank):
        """
        Calculates the rank of End(G) based on the extension type, using the formula:
        rank(End(G)) = rank(Hom(G, T)) + rank(End(A,[G]))
        where End(A,[G]) is the subring of End(A) preserving the extension class.
        """
        rank_End_A = self.A.endomorphism_ring.rank
        rank_End_T = self.T.endomorphism_ring.rank

        if self.extension_type == 'split':
            # For a split extension G = A x T:
            # End(A,[G]) = End(A) because the extension class is trivial.
            rank_End_A_preserving_extension = rank_End_A
            # Hom(G,T) = Hom(A x T, T) is isomorphic to Hom(T, T) = End(T).
            rank_Hom_G_T = rank_End_T
            rank_End_G = rank_Hom_G_T + rank_End_A_preserving_extension
            equation = f"rank(End(T)) + rank(End(A)) = {rank_Hom_G_T} + {rank_End_A_preserving_extension} = {rank_End_G}"
            return EndomorphismRing(rank_End_G), equation

        elif self.extension_type == 'rigid_non_split':
            # For a "generic" or "rigid" non-split extension, we can model a plausible scenario:
            # 1. End(A,[G]) is small. Model its rank as rigid_rank (e.g., 1, for just Z).
            rank_End_A_preserving_extension = rigid_rank
            # 2. Hom(G,T) is small. Its rank is at most rank(End(T)).
            #    Let's assume its rank is also small, e.g. min(rigid_rank, rank_End_T)
            rank_Hom_G_T = min(rigid_rank, rank_End_T)
            rank_End_G = rank_Hom_G_T + rank_End_A_preserving_extension
            equation = f"rank(Hom(G,T)) + rank(End(A,[G])) ≈ {rank_Hom_G_T} + {rank_End_A_preserving_extension} = {rank_End_G}"
            return EndomorphismRing(rank_End_G), equation
        else:
            raise ValueError("Unknown extension type")

# ---- Illustrate two scenarios to demonstrate the ambiguity ----

# Scenario 1: Split case, where G has more endomorphisms.
print("--- Scenario 1: Split Extension (G = A x T) ---")
A1 = AbelianVariety(name="A1 (typical elliptic curve)", endomorphism_rank=1)
T1 = Torus(dimension=1)
G1 = SemiAbelianVariety(A1, T1, extension_type='split')
print(f"Abelian Variety A1 has End(A1) rank = {A1.endomorphism_ring.rank}.")
print(f"Semi-Abelian Variety G1 (a split extension) has End(G1) rank calculated as:")
print(f"rank(End(G1)) = {G1.equation}")
if G1.endomorphism_ring.rank > A1.endomorphism_ring.rank:
    print("Result: In this scenario, G has more endomorphisms than A.")

# Scenario 2: Rigid non-split case, where A has more endomorphisms.
print("\n--- Scenario 2: Rigid Non-Split Extension ---")
A2 = AbelianVariety(name="A2 (CM abelian surface)", endomorphism_rank=8)
T2 = Torus(dimension=1)
G2 = SemiAbelianVariety(A2, T2, extension_type='rigid_non_split', rigid_rank=1)
print(f"Abelian Variety A2 has End(A2) rank = {A2.endomorphism_ring.rank}.")
print(f"Semi-Abelian Variety G2 (a rigid extension) has End(G2) rank calculated based on a model:")
print(f"rank(End(G2)) = {G2.equation}")
if A2.endomorphism_ring.rank > G2.endomorphism_ring.rank:
    print("Result: In this scenario, A has more endomorphisms than G.")

print("\n--- Conclusion ---")
print("As the two scenarios show, the relative number of endomorphisms depends on the specific extension.")
