def solve():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.
    This is equivalent to enumerating the dual graphs of stable curves of
    genus 3 with 1 marked point and 2 nodes.
    """
    print("This program calculates the number of codimension 2 boundary strata of the moduli space of stable genus 3 curves with 1 marked point, M_bar_{3,1}.")
    print("This is equivalent to counting the number of non-isomorphic stable dual graphs with 2 edges, representing curves of arithmetic genus 3 with 1 marked point.\n")

    def is_stable(g, k):
        """Checks the stability condition 2g - 2 + k > 0."""
        return 2 * g - 2 + k > 0

    total_strata = 0
    strata_counts = []

    print("--- Case 1: Irreducible Curves (1 vertex, 2 loops) ---")
    # |V|=1, |E|=2 => h1 = 2-1+1 = 2.
    # Genus equation: g_total = g_v + h1 => 3 = g_v + 2 => g_v = 1.
    # The single component has g=1. It has 2 nodes and 1 marked point, so k=3.
    # Stability: 2*1 - 2 + 3 = 3 > 0. This is stable.
    # This configuration is unique.
    count1 = 1
    total_strata += count1
    strata_counts.append(count1)
    print("A single component of genus 1 with 2 nodes and 1 marked point.")
    print(f"Number of strata in this case: {count1}\n")


    print("--- Case 2: Annulus-type Curves (2 vertices, 2 edges forming a loop) ---")
    # |V|=2, |E|=2 => h1 = 2-2+1 = 1.
    # Genus equation: g1 + g2 + h1 = 3 => g1 + g2 = 2.
    count2 = 0
    # Subcase 2.1: Genera (2, 0).
    # To be stable, the g=0 component must have at least 3 special points.
    # The two nodes provide k=2. The marked point must be on the g=0 component, making k=3.
    # v1(g=2, k=2 nodes) is stable. v2(g=0, k=2 nodes + 1 mp = 3) is stable.
    count2 += 1
    print("A component of genus 2 and a component of genus 0, meeting at two points. The marked point is on the genus 0 component.")
    # Subcase 2.2: Genera (1, 1).
    # Both components have genus 1. They are symmetric. The marked point can be on either.
    # v1(g=1, k=2 nodes + 1 mp = 3) is stable. v2(g=1, k=2 nodes) is stable.
    count2 += 1
    print("Two components of genus 1, meeting at two points. The marked point is on one of them.")
    total_strata += count2
    strata_counts.append(count2)
    print(f"Number of strata in this case: {count2}\n")


    print("--- Case 3: Lollipop-type Curves (2 vertices, 1 edge, 1 loop) ---")
    # |V|=2, |E|=2 => h1 = 1. Genus equation: g1 + g2 = 2.
    # One vertex (v1) has a self-loop and an edge to v2. k1 >= 2, k2 >= 1.
    count3 = 0
    # Subcase 3.1: Genera (1, 1). Let v1 have the loop.
    # MP on v1: v1(g=1, k=3), v2(g=1, k=1). Both stable.
    count3 += 1
    print("Two components of genus 1. One has a self-node, and the marked point is on it. The other is attached to it.")
    # MP on v2: v1(g=1, k=2), v2(g=1, k=2). Both stable.
    count3 += 1
    print("Two components of genus 1. One has a self-node, and the marked point is on the other attached component.")
    # Subcase 3.2: Genera (2, 0).
    # Loop on g=0 component (v1): MP must be on v1 for stability (k=3). v1(g=0, k=3), v2(g=2, k=1). Stable.
    count3 += 1
    print("A component of genus 0 with a self-node and the marked point, attached to a component of genus 2.")
    # Other combinations with (g=0, g=2) are unstable (e.g., loop on g=0, mp on g=2 leads to k=2 for g=0 component).
    total_strata += count3
    strata_counts.append(count3)
    print(f"Number of strata in this case: {count3}\n")


    print("--- Case 4: Chain-type Curves (3 components in a line) ---")
    # |V|=3, |E|=2 => h1 = 0.
    # Genus equation: g1 + g2 + g3 = 3.
    count4 = 0
    # Subcase 4.1: Marked point on the central component (v2).
    # k1=1, k2=3, k3=1. Stability requires g1>=1, g3>=1.
    #   - Genera (1,1,1): v1(g=1,k=1), v2(g=1,k=3), v3(g=1,k=1). Stable.
    count4 += 1
    print("A chain of three genus 1 components, with the marked point on the central one.")
    #   - Genera (1,0,2) or (2,0,1): v1(g=1,k=1), v2(g=0,k=3), v3(g=2,k=1). Stable.
    count4 += 1
    print("A chain of components of genera (1, 0, 2), with the marked point on the central (genus 0) one.")
    # Subcase 4.2: Marked point on an end component (v1).
    # k1=2, k2=2, k3=1. Stability requires g1>=1, g2>=1, g3>=1.
    # The only partition of 3 into 3 positive integers is (1,1,1).
    count4 += 1
    print("A chain of three genus 1 components, with the marked point on an end one.")
    total_strata += count4
    strata_counts.append(count4)
    print(f"Number of strata in this case: {count4}\n")


    print("--- Total Calculation ---")
    equation = " + ".join(map(str, strata_counts))
    print(f"The total number of strata is the sum of the counts from each case:")
    print(f"{equation} = {total_strata}")

solve()
<<<9>>>