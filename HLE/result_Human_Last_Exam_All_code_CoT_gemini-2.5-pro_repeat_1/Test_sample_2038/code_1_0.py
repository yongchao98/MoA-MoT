# The user needs to have the 'snappy' and 'sympy' libraries installed.
# They can be installed using pip:
# pip install snappy sympy

import snappy
from sympy.ntheory import sqrt_mod, primefactors

def find_knots():
    """
    Finds and counts 2-bridge knots with crossing number at most 13
    that admit two disjoint non-parallel embedded minimal genus Seifert surfaces.
    """
    # The problem specifies a parameter 'a' >= 2, which corresponds to 'n'
    # in the relevant Alexander polynomial formula.
    n_min = 2

    # We need to decide how far to check for n. The crossing number c(K(p,q))
    # tends to grow with p = 4n+1. A few hundred should be sufficient to find
    # all knots with a crossing number up to 13.
    n_max = 300

    found_knots = set()
    equation_terms = []

    # Loop through possible values of n
    for n in range(n_min, n_max + 1):
        p = 4 * n + 1

        # Check if solutions for q^2 = -1 (mod p) exist.
        # This is true if and only if all prime factors of p are of the form 4k+1.
        if any(factor % 4 == 3 for factor in primefactors(p)):
            continue

        # Find all solutions for q. sympy.sqrt_mod handles the number theory.
        solutions = sqrt_mod(-1, p, all_roots=True)

        if not solutions:
            continue

        # To count a knot and its mirror image once, we take one representative 'q'
        # from each pair {q, p-q}. We choose the smaller one.
        q_candidates = sorted(list({min(q, p - q) for q in solutions}))

        for q in q_candidates:
            if q == 0:
                continue

            try:
                # Create the 2-bridge knot object
                knot = snappy.Link(twobridge=[p, q])

                # We are looking for knots (1-component links)
                if knot.num_components() == 1:
                    # Check if we have already found this knot
                    if knot in found_knots:
                        continue

                    # Calculate the crossing number
                    cn = knot.crossing_number()

                    if cn <= 13:
                        found_knots.add(knot)
                        equation_terms.append(1)
                        print(f"Found a knot from n={n} (p={p}, q={q}) with crossing number {cn}. This is knot {knot.identify()}.")

            except Exception:
                # This may happen for invalid (p,q) pairs, though our method should avoid them.
                # We can safely ignore these cases.
                pass

    print("\nCalculating the total number of such knots:")
    equation_str = " + ".join(map(str, equation_terms))
    print(f"{equation_str} = {len(found_knots)}")
    
    return len(found_knots)

if __name__ == '__main__':
    final_count = find_knots()
    print(f"\nIn total, there are {final_count} such knots.")
    print(f"<<<{final_count}>>>")
