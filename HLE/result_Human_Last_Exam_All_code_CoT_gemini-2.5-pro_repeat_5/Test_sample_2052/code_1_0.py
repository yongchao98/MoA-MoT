import math

def calculate_s_lambda_at_delta_kn(lam, n, char_table, class_indices, z_values):
    """
    Calculates s_lambda(delta_k,n) which is chi_lambda(mu) / z_mu,
    where mu is the partition of |lambda| into n's.
    """
    size = sum(lam)
    if size % n != 0:
        return 0
    
    num_cycles = size // n
    class_partition_str = str(n)
    
    # This simplified example only handles cases where the cycle type exists in the pre-computed tables.
    if class_partition_str not in class_indices:
        # For this example, if the required cycle type for n is not in the S4 table, we can't compute it.
        # e.g., for n=3, we would need cycle type (3,1), not (3,3,...)
        # The logic holds for any |lambda|, but we demonstrate with |lambda|=4.
        return float('nan')

    class_idx = class_indices[class_partition_str]
    chi_val = char_table[lam][class_idx]
    z_val = z_values[class_partition_str]
    
    return chi_val / z_val


def main():
    """
    Demonstrates that Z_-4 is not trivially related to Z_-2.
    """
    print("--- Analysis for Question (b) ---")
    print("We check if the ratio of coefficients s_lambda(delta_k,4) / s_lambda(delta_k,2) is constant.")
    print("The check is performed for all partitions lambda of size 4.\n")

    # Partitions of 4
    partitions_of_4 = [
        (4,),
        (3, 1),
        (2, 2),
        (2, 1, 1),
        (1, 1, 1, 1),
    ]

    # Character table for the Symmetric Group S4
    # Source: Fulton and Harris, "Representation Theory: A First Course", p. 19
    # Rows are irreps (partitions), columns are conjugacy classes.
    # Conjugacy classes order: (1,1,1,1), (2,1,1), (3,1), (4), (2,2)
    char_table_S4 = {
        (4,):           [1, 1, 1, 1, 1],
        (3, 1):         [3, 1, 0, -1, -1],
        (2, 2):         [2, 0, -1, 0, 2],
        (2, 1, 1):      [3, -1, 0, 1, -1],
        (1, 1, 1, 1):   [1, -1, 1, -1, 1],
    }
    
    # To compute s_lambda(delta_k,n), we need characters for cycle types (n,n,...)
    # For |lambda|=4, the relevant cycle types are (4) and (2,2).
    class_indices_S4 = {'4': 3, '2,2': 4}

    # z_mu = product(k^m_k * m_k!) where mu = 1^m_1 2^m_2 ...
    z_mu_values = {
        '4': 4,      # for class (4), z = 4^1 * 1! = 4
        '2,2': 8,    # for class (2,2), z = 2^2 * 2! = 8
    }

    print(f"{'Partition':<12} | {'s_lam(d_k,2)':<12} | {'s_lam(d_k,4)':<12} | {'Ratio':<12}")
    print("-" * 55)

    for p in partitions_of_4:
        # For n=2, |lambda|=4, cycle type is (2,2)
        s_lam_d2 = calculate_s_lambda_at_delta_kn(p, 2, char_table_S4, {'2,2': 0}, {'2,2': 8})
        # For this simplified implementation, we need to adjust inputs for the generic function
        chi_val_d2 = char_table_S4[p][class_indices_S4['2,2']]
        s_lam_d2 = chi_val_d2 / z_mu_values['2,2']
        
        # For n=4, |lambda|=4, cycle type is (4)
        chi_val_d4 = char_table_S4[p][class_indices_S4['4']]
        s_lam_d4 = chi_val_d4 / z_mu_values['4']

        ratio = "undefined"
        if s_lam_d2 != 0:
            ratio = s_lam_d4 / s_lam_d2
            ratio = f"{ratio:.2f}"

        print(f"{str(p):<12} | {s_lam_d2:<12.3f} | {s_lam_d4:<12.3f} | {ratio:<12}")

    print("\nConclusion: The ratio is not constant. Therefore, Z_-4 and Z_-2 are not")
    print("trivially proportional, and a specific operator W_-4 is necessary to generate Z_-4.")

if __name__ == "__main__":
    main()