import math

def solve():
    """
    Calculates the number of Tverberg partitions for a special point set configuration
    based on a given integer partition.
    """
    # The user should modify this line to provide the partition a_1, ..., a_r.
    # Example partition:
    partition_parts = [3, 2, 1]

    r = len(partition_parts)
    n = sum(partition_parts)

    print(f"The given partition of n={n} is a = {partition_parts}.")
    print(f"The number of parts is r = {r}.")
    
    # Formula for the number of labeled Tverberg partitions:
    # N = (r-1)! * (a_1-1)! * (a_2-1)! * ... * (a_r-1)!
    
    # Calculate (r-1)!
    try:
        r_minus_1_fact = math.factorial(r - 1)
    except ValueError:
        # This case happens if r=0, i.e., empty partition list.
        r_minus_1_fact = 1

    # Calculate product of (a_i - 1)!
    prod_a_i_minus_1_fact = 1
    for part_size in partition_parts:
        if part_size < 1:
            print("Error: All part sizes must be at least 1.")
            return
        prod_a_i_minus_1_fact *= math.factorial(part_size - 1)

    # Calculate the total number of partitions
    num_partitions = r_minus_1_fact * prod_a_i_minus_1_fact

    # Print the step-by-step calculation
    print("\nThe number of Tverberg partitions is calculated using the formula:")
    print("N = (r-1)! * (a_1-1)! * (a_2-1)! * ... * (a_r-1)!")
    
    # Build the string for the formula with numbers
    formula_str = f"N = ({r}-1)! * "
    part_calcs = []
    for part_size in partition_parts:
        part_calcs.append(f"({part_size}-1)!")
    formula_str += " * ".join(part_calcs)
    print(formula_str)

    # Build the string for the factorial values
    fact_val_str = f"N = {r_minus_1_fact} * "
    part_fact_vals = []
    for part_size in partition_parts:
        part_fact_vals.append(str(math.factorial(part_size - 1)))
    fact_val_str += " * ".join(part_fact_vals)
    print(fact_val_str)
    
    # Print the final result
    print(f"N = {num_partitions}")

    print(f"\nFinal Answer: {num_partitions}")

solve()