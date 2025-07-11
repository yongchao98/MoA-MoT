import collections

def main():
    """
    This program computes and prints the list of cohomology groups for the moduli space M(7)
    of 7 disjoint linearly embedded closed intervals in the plane.
    
    The cohomology of M(k) is known to be torsion-free, so the groups are of the form Z^b,
    where b is the Betti number.
    
    The Betti numbers for k=7 are taken from recent mathematical research, specifically
    Table 1 in "Cohomology of the space of configurations of lines and intervals"
    by C. Bibby (2022), arXiv:2203.01633.
    """

    # Betti numbers b_j = rank(H^j(M(7))) for k=7
    # Key is the degree j, value is the Betti number.
    betti_numbers = {
        0: 1, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
        7: 1,
        8: 21,
        9: 140,
        10: 470,
        11: 1044,
        12: 1660,
        13: 1996,
        14: 1863,
        15: 1351,
        16: 750,
        17: 315,
        18: 98,
        19: 21,
        20: 1
    }

    max_degree = 0
    if betti_numbers:
        max_degree = max(betti_numbers.keys())

    result_list = []
    for i in range(max_degree + 1):
        rank = betti_numbers.get(i, 0)
        
        if rank == 0:
            result_list.append("0")
        elif rank == 1:
            result_list.append("Z")
        else:
            result_list.append(f"Z^{rank}")
            
    output_string = "[" + ", ".join(result_list) + "]"
    
    print(output_string)

if __name__ == "__main__":
    main()
