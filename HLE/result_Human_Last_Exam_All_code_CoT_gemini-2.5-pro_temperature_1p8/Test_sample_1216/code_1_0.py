import math

def solve_augustus_records():
    """
    Solves the Augustus records problem by finding the integer solution to a system of equations
    derived from the problem's constraints.
    """
    solution = None
    # We determined that d = k*k, and f, s are dependent on k and l.
    # We derived the equation: 16k^2 + 32k + l^2 + 96l - 9600 = 0
    # Let's iterate through possible values to find the integer solution.
    # d must be smaller than readable records (600), so k < sqrt(600) approx 24.
    for k in range(1, 25):
        # l must be a multiple of 4. l must be less than 120 and also s = (l/4)^2 < 600 -> l < 98.
        for l in range(4, 97, 4):
            # Check if the values satisfy the equation
            if 16*k**2 + 32*k + l**2 + 96*l - 9600 == 0:
                solution = (k, l)
                break
        if solution:
            break

    if not solution:
        print("No integer solution found within the constraints.")
        return

    k, l = solution
    
    # Calculate the four key quantities based on the solution
    # a) Lost records (Octavius ones)
    lost_records_octavius = l
    
    # b) Documents with dual naming (Augustus+Caesar)
    dual_named_records = k**2
    
    # c) Single-variant documents
    single_variant_records = (l // 4)**2
    
    # d) Full imperial title documents
    full_title_records = k + 3 * l

    # --- Verification Step ---
    # The sum of all numbers found must be divisible by the number of distinct naming patterns.
    caesar_only = 80
    total_readable = 600
    
    # The five mutually exclusive categories of readable records are:
    # dual_named, single_variant, full_title, caesar_only, and icdf.
    # We can find the size of the 'icdf' category by subtraction.
    icdf_records = total_readable - (dual_named_records + single_variant_records + full_title_records + caesar_only)
    
    # Sum of all key numbers found (the 5 readable categories + the octavius-lost ones)
    sum_of_numbers = dual_named_records + single_variant_records + full_title_records + caesar_only + icdf_records + lost_records_octavius
    
    # The number of distinct naming patterns identified in the readable archive
    num_distinct_patterns = 5
    
    if sum_of_numbers % num_distinct_patterns != 0:
        print(f"Verification failed: Sum of numbers ({sum_of_numbers}) is not divisible by the number of patterns ({num_distinct_patterns}).")
        return
        
    # --- Final Calculation ---
    
    # The product of the four main calculated values
    product = lost_records_octavius * dual_named_records * single_variant_records * full_title_records
    
    # Denominator is "The sum of the distinct ways Augustus is named", which we interpret as the number of distinct patterns
    denominator = num_distinct_patterns
    
    answer = product / denominator
    final_answer = math.floor(answer / 1000)
    
    print("Final Equation Components:")
    print(f"a) Lost records (Octavius): {lost_records_octavius}")
    print(f"b) Dual-named records: {dual_named_records}")
    print(f"c) Single-variant records: {single_variant_records}")
    print(f"d) Full imperial title records: {full_title_records}")
    print(f"Denominator (Distinct naming patterns): {denominator}")
    
    print("\nCalculation:")
    print(f"({lost_records_octavius} * {dual_named_records} * {single_variant_records} * {full_title_records}) / {denominator} = {answer}")
    
    print("\nFinal Answer (floor(answer / 1000)):")
    print(f"floor({answer} / 1000) = {final_answer}")
    
    # Per instruction, output the final answer in the specified format
    print(f"<<<{final_answer}>>>")


solve_augustus_records()