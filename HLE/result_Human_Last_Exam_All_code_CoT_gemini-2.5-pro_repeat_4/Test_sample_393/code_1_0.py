import math
from collections import Counter

def solve_tverberg_partitions():
    """
    Calculates the number of Tverberg partitions for a given integer partition.
    
    The problem asks for the number of ways to partition a set of n points,
    given that all partitions must have subsets of sizes corresponding to a
    given integer partition a = (a_1, ..., a_r).
    
    The formula is: (n! / (a_1! * a_2! * ... * a_r!)) / (k_1! * k_2! * ...),
    where k_i is the number of times a particular part size appears.
    """
    
    # Example partition a = (a_1, a_2, ..., a_r)
    a = [4, 3, 3, 2]
    
    # Step 1: Calculate n, the total number of points.
    n = sum(a)
    
    print(f"Given the integer partition a = {a}")
    print("-" * 30)
    
    # Print the calculation of n
    n_calculation_str = " + ".join(map(str, a))
    print(f"Step 1: Calculate the total number of points, n.")
    print(f"n = {n_calculation_str} = {n}\n")
    
    # Step 2: Calculate the multinomial coefficient for ordered partitions.
    # M = n! / (a_1! * a_2! * ... * a_r!)
    print("Step 2: Calculate the number of ordered partitions (multinomial coefficient M).")
    print(f"M = n! / ({' * '.join([f'{i}!' for i in a])})")
    
    numerator = math.factorial(n)
    denominator_parts = [math.factorial(i) for i in a]
    denominator = 1
    for part in denominator_parts:
        denominator *= part
        
    print(f"M = {n}! / ({' * '.join(map(str, denominator_parts))})")
    print(f"M = {numerator} / {denominator}")
    
    m_result = numerator // denominator
    print(f"M = {m_result}\n")

    # Step 3: Count frequencies of each part size to correct for unordered partitions.
    print("Step 3: Correct for unordered partitions by counting part size frequencies.")
    size_counts = Counter(a)
    print(f"The frequencies of the part sizes are: {dict(size_counts)}")
    
    # Step 4: Calculate the correction factor D.
    # D = k_1! * k_2! * ...
    print("Step 4: Calculate the correction factor D from these frequencies.")
    correction_denominator_parts = [math.factorial(count) for count in size_counts.values()]
    correction_denominator = 1
    for part in correction_denominator_parts:
        correction_denominator *= part
    
    correction_calculation_str = ' * '.join([f'{count}!' for count in size_counts.values()])
    correction_values_str = ' * '.join(map(str, correction_denominator_parts))
    print(f"D = {correction_calculation_str}")
    if len(correction_denominator_parts) > 1:
        print(f"D = {correction_values_str} = {correction_denominator}\n")
    else:
        print(f"D = {correction_denominator}\n")

    # Step 5: Final calculation.
    # Total = M / D
    final_result = m_result // correction_denominator
    print("Step 5: Calculate the final number of Tverberg partitions.")
    print(f"Total Partitions = M / D")
    print(f"Total Partitions = {m_result} / {correction_denominator}")
    print(f"Total Partitions = {final_result}")
    
    return final_result

# Run the solver and capture the final answer
final_answer = solve_tverberg_partitions()
print(f"\n<<< {final_answer} >>>")
