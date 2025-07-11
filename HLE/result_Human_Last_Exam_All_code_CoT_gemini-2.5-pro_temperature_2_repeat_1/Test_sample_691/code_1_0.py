import sys

# On some platforms, the 'x' character may not render correctly in the console.
# Using a simple try-except block to handle encoding for wider compatibility.
try:
    sys.stdout.reconfigure(encoding='utf-8')
except TypeError:
    pass # Not all environments support this (e.g., in-browser notebooks)

def solve_topology_problem():
    """
    This function provides a step-by-step derivation for the fundamental group
    of the described topological space.
    """
    print("### Step 1: Identify the initial surface M ###")
    print("A pair of pants is topologically a sphere with 3 holes (a 3-holed sphere), or S_{g,b} where g=0, b=3.")
    
    g_pants = 0
    b_pants = 3
    # The Euler characteristic chi is given by chi = 2 - 2g - b
    chi_pants = 2 - 2 * g_pants - b_pants
    print(f"The Euler characteristic for one pair of pants is: 2 - 2*{g_pants} - {b_pants} = {chi_pants}")

    # We sew two pairs of pants (P1, P2) along two leg openings each. The seam consists of two circles.
    # The Euler characteristic of a circle is 0.
    # chi(M) = chi(P1) + chi(P2) - chi(seam)
    chi_M = chi_pants + chi_pants - 0
    print(f"The surface M, made by sewing two pants, has Euler characteristic: {chi_pants} + {chi_pants} - 0 = {chi_M}")

    # The resulting surface M has 2 boundary components (the two waistbands).
    b_M = 2
    # We find the genus g_M of M using the formula again: chi = 2 - 2g - b
    # -2 = 2 - 2*g_M - 2  =>  -2 = -2*g_M  => g_M = 1
    g_M = (2 - b_M - chi_M) / 2
    print(f"M has {b_M} boundaries. Its genus g can be found from {chi_M} = 2 - 2*g - {b_M}, which means g = {int(g_M)}.")
    print("Thus, M is a genus-1 surface with 2 boundaries (a torus with two holes), S_{1,2}.")
    print("-" * 40)

    print("### Step 2: Find the fundamental group of M ###")
    # The fundamental group of a surface S_{g,b} with b>0 is the free group F_n on n generators.
    # The rank n is given by n = 2g + b - 1
    rank = 2 * g_M + b_M - 1
    print(f"The rank of the fundamental group for S_{{{int(g_M)}, {b_M}}} is n = 2*g + b - 1.")
    print(f"The calculation is: n = 2*{int(g_M)} + {b_M} - 1 = {int(rank)}")
    print(f"So, the fundamental group of M is the free group on {int(rank)} generators, F_{int(rank)}.")
    print("-" * 40)

    print("### Step 3: Add relations from collapsing the waistbands ###")
    print(f"Let the fundamental group of M be <a, b, c>, corresponding to F_{int(rank)}.")
    print("Collapsing the two boundary loops (the waistbands) to a point makes them trivial in the group.")
    print("The boundary loops can be expressed in terms of the generators as c and [a,b]c (where [a,b] = aba⁻¹b⁻¹).")
    print("This imposes two relations:")
    print("1. c = 1")
    print("2. [a,b]c = 1")
    print("Substituting the first relation into the second gives [a,b] = 1.")
    print("The relation [a,b]=1 is equivalent to the commutator relation ab = ba.")
    print("-" * 40)

    print("### Step 4: Final group identification ###")
    print("The final group has the presentation <a, b | ab=ba>.")
    print("This is the free abelian group on two generators, which is the direct product of Z with itself.")
    
    # Final equation format as requested
    print("\nFinal Fundamental Group Equation:")
    final_eq_parts = ["Z", "x", "Z"]
    print(" ".join(final_eq_parts))

solve_topology_problem()
