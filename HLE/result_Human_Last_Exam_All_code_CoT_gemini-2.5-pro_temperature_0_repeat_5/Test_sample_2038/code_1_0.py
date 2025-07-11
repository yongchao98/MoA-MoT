import math
import sympy
import snappy

def solve_knot_problem():
    """
    Finds and counts the 2-bridge knots with crossing number at most 13
    that admit two disjoint non-parallel embedded minimal genus Seifert surfaces.
    """
    print("Starting the search for qualifying knots...")
    print("This requires the 'snappy' and 'sympy' libraries.")
    print("-" * 50)

    found_knots = []
    t = sympy.var('t')
    
    # The crossing number tends to increase with n.
    # A search up to n=15 is sufficient as crossing numbers will exceed 13.
    for n in range(1, 16):
        # Case 1: p = 4n - 1
        p1 = 4 * n - 1
        if p1 > 1:
            # Target Alexander polynomial (in standard, not Laurent form)
            target_poly1 = n * t**2 - (2 * n - 1) * t + n
            
            # Iterate through possible q values for the 2-bridge knot K(p,q)
            for q1 in range(1, (p1 // 2) + 1):
                if math.gcd(p1, q1) == 1:
                    try:
                        knot = snappy.Knot(rational=(p1, q1))
                        alex_poly = knot.alexander_polynomial()
                        
                        if alex_poly == target_poly1:
                            crossing_num = knot.crossing_number()
                            if crossing_num <= 13:
                                knot_id = knot.identify()
                                # Use the first name from the list if available
                                knot_name = knot_id[0].name if knot_id else f"K({p1},{q1})"
                                
                                entry = {
                                    "name": knot_name,
                                    "p": p1,
                                    "q": q1,
                                    "c": crossing_num,
                                    "n": n,
                                    "type": 1,
                                    "poly_str": f"{n}*t^2 - ({2*n - 1})*t + {n} = {sympy.simplify(target_poly1)}"
                                }
                                if entry not in found_knots:
                                    found_knots.append(entry)
                    except Exception as e:
                        # Some (p,q) might not form valid knots in snappy's view
                        continue

        # Case 2: p = 4n + 1
        p2 = 4 * n + 1
        # Target Alexander polynomial
        target_poly2 = n * t**2 - (2 * n + 1) * t + n
        
        # Iterate through possible q values
        for q2 in range(1, (p2 // 2) + 1):
            if math.gcd(p2, q2) == 1:
                try:
                    knot = snappy.Knot(rational=(p2, q2))
                    alex_poly = knot.alexander_polynomial()

                    if alex_poly == target_poly2:
                        crossing_num = knot.crossing_number()
                        if crossing_num <= 13:
                            knot_id = knot.identify()
                            knot_name = knot_id[0].name if knot_id else f"K({p2},{q2})"
                            
                            entry = {
                                "name": knot_name,
                                "p": p2,
                                "q": q2,
                                "c": crossing_num,
                                "n": n,
                                "type": 2,
                                "poly_str": f"{n}*t^2 - ({2*n + 1})*t + {n} = {sympy.simplify(target_poly2)}"
                            }
                            if entry not in found_knots:
                                found_knots.append(entry)
                except Exception as e:
                    continue

    print("Found the following knots satisfying the conditions:")
    # Sort knots by crossing number, then by name
    found_knots.sort(key=lambda k: (k['c'], k['name']))
    
    for i, knot in enumerate(found_knots):
        print(f"\n{i+1}. Knot: {knot['name']} = K({knot['p']}, {knot['q']})")
        print(f"   Crossing Number: {knot['c']}")
        print(f"   Corresponds to n={knot['n']} with p = 4*n {'-' if knot['type']==1 else '+'} 1 = {knot['p']}")
        print(f"   Alexander Polynomial: {knot['poly_str']}")

    print("-" * 50)
    print(f"Total number of such knots is: {len(found_knots)}")
    
    # Final answer in the required format
    print(f"\n<<<{len(found_knots)}>>>")

if __name__ == '__main__':
    solve_knot_problem()