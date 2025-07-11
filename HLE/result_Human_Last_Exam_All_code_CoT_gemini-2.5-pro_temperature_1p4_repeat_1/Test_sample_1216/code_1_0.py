import math

def solve_augustus_records():
    """
    Solves the Augustus records problem by calculating each variable step-by-step based on the prompt's constraints.
    """
    # Step 1: Initial Data Extraction
    total_records = 720
    total_lost_records = total_records / 6
    readable_records = total_records - total_lost_records
    
    # These are fixed categories within the readable records
    icdf_records = readable_records / 3
    caesar_only_records = 80
    
    # Step 2: Use constraints to find the core variables
    
    # The prompt differentiates between "the number of lost records" (total) and "the Octavius-lost records" (a specific part).
    # This is crucial for solving the constraints.
    
    # "The records using single variants... root equals the number of lost records divided by 4"
    # This refers to the total number of lost records.
    s_sqrt = total_lost_records / 4
    s_single_variant_docs = s_sqrt**2
    
    # Now we establish a relationship between D, F, and L (Octavius-lost).
    # The sum of all readable categories must equal the total readable records.
    # The given categories are I, C, D, F. We assume these are the only ones.
    # D + F + icdf_records + caesar_only_records = readable_records
    # D + F = readable_records - icdf_records - caesar_only_records
    d_plus_f_sum = readable_records - icdf_records - caesar_only_records
    
    # We now have:
    # 1) F = sqrt(D) + 3*L
    # 2) D + F = d_plus_f_sum
    # Substituting F from (1) into (2):
    # D + sqrt(D) + 3*L = d_plus_f_sum
    # We need to find L (Octavius-lost) and D (Dual-named)
    # We know L must be a multiple of 4 (from another constraint in the full text logic),
    # and D must be a perfect square. We can iterate to find a solution.
    
    l_octavius_lost = 0
    d_dual_named_docs = 0
    
    # Iterate through possible values for L (must be a multiple of 4 as per intermediate logic not shown here for brevity)
    # L is a part of total_lost_records, so L <= 120
    for l_candidate in range(4, int(total_lost_records) + 1, 4):
        # We need to solve for D in: D + sqrt(D) = d_plus_f_sum - 3*L
        target = d_plus_f_sum - 3 * l_candidate
        
        # Solving the quadratic equation x^2 + x - target = 0 for x = sqrt(D)
        discriminant = 1 + 4 * target
        if discriminant >= 0 and math.isqrt(discriminant)**2 == discriminant:
            # Check if it gives an integer solution for sqrt(D)
            sqrt_d = (-1 + math.isqrt(discriminant)) / 2
            if sqrt_d == int(sqrt_d) and sqrt_d > 0:
                l_octavius_lost = l_candidate
                d_dual_named_docs = int(sqrt_d**2)
                break
                
    # Now calculate F using the found values
    f_full_title_docs = math.sqrt(d_dual_named_docs) + 3 * l_octavius_lost
    
    # Step 3: Identify the number of distinct naming patterns for "Augustus" in the archive
    # The archive consists of readable records. S=900 is not a readable record, so it's not in the archive.
    # 1. "Augustus" and "Caesar" (D) -> Yes, contains "Augustus"
    # 2. "Imperator Caesar Augustus" (F) -> Yes, contains "Augustus"
    # The other readable types (ICDF, C) do not contain the name "Augustus".
    distinct_patterns = 2

    # Step 4: Perform the final calculation
    # Product of: a) Lost records (Octavius) b) Dual naming c) Single-variant d) Full imperial title
    product_of_values = l_octavius_lost * d_dual_named_docs * s_single_variant_docs * f_full_title_docs

    # Divided by: The sum of the distinct ways Augustus is named
    answer = product_of_values / distinct_patterns
    
    final_answer = math.floor(answer / 1000)

    # Output the steps and the final answer
    print("Step 1: Calculating the values of the four key document types.")
    print(f"a) The number of Octavius-lost records (L) is: {int(l_octavius_lost)}")
    print(f"b) The number of dual-named documents (D) is: {int(d_dual_named_docs)}")
    print(f"c) The number of single-variant documents (S) is: {int(s_single_variant_docs)}")
    print(f"d) The number of full imperial title documents (F) is: {int(f_full_title_docs)}")
    print("\nStep 2: Performing the final calculation.")
    print(f"The calculation is: (a * b * c * d) / (Number of Distinct Naming Patterns)")
    print(f"Which translates to: ({int(l_octavius_lost)} * {int(d_dual_named_docs)} * {int(s_single_variant_docs)} * {int(f_full_title_docs)}) / {distinct_patterns}")
    print(f"Result of the calculation = {answer}")
    print(f"\nStep 3: The final answer is floor(Result / 1000).")
    print(f"Final Answer = floor({answer} / 1000) = {final_answer}")
    
    return final_answer

final_result = solve_augustus_records()
print(f"<<<{final_result}>>>")