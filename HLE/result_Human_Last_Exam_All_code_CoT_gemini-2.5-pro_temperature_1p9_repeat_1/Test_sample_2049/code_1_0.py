import math

def main():
    """
    Calculates the number of isomorphism classes of del Pezzo surfaces of degree 5
    over Q with good reduction everywhere except possibly at the prime 2.
    """
    
    print("This problem is equivalent to counting the number of étale Q-algebras of dimension 4 that are unramified outside the prime 2.")
    print("We do this by counting the number of ways to form such an algebra from products of number fields that are unramified outside 2.\n")

    # Step 1: Define the number of building blocks (number fields unramified outside 2)
    # k_n is the number of number fields of degree n unramified outside 2.
    k1 = 1  # The field Q
    k2 = 3  # The fields Q(i), Q(sqrt(2)), Q(sqrt(-2))
    k3 = 0  # No such cubic fields exist
    k4 = 7  # The number of such quartic fields is a known result from number theory

    print("--- Number of Building Block Fields ---")
    print(f"Number of degree 1 fields unramified outside 2: {k1}")
    print(f"Number of degree 2 fields unramified outside 2: {k2}")
    print(f"Number of degree 3 fields unramified outside 2: {k3}")
    print(f"Number of degree 4 fields unramified outside 2: {k4}\n")

    print("--- Counting Algebras by Partition of 4 ---")
    
    # Partition [4]: One field of degree 4
    num_partition_4 = k4
    print(f"Type [4] (one quartic field): {num_partition_4}")

    # Partition [3,1]: One field of degree 3, one of degree 1
    num_partition_3_1 = k3 * k1
    print(f"Type [3,1] (cubic x linear): {k3} * {k1} = {num_partition_3_1}")
    
    # Partition [2,2]: Two fields of degree 2
    # Case A: Two distinct quadratic fields, C(k2, 2)
    num_p22_distinct = math.comb(k2, 2)
    # Case B: One quadratic field squared, C(k2, 1)
    num_p22_same = math.comb(k2, 1)
    num_partition_2_2 = num_p22_distinct + num_p22_same
    print(f"Type [2,2] (quadratic x quadratic): {num_p22_distinct} (distinct) + {num_p22_same} (same) = {num_partition_2_2}")

    # Partition [2,1,1]: One field of degree 2, two of degree 1
    num_partition_2_1_1 = k2
    print(f"Type [2,1,1] (quadratic x linear x linear): {num_partition_2_1_1}")

    # Partition [1,1,1,1]: Four fields of degree 1
    num_partition_1_1_1_1 = 1
    print(f"Type [1,1,1,1] (linear x linear x linear x linear): {num_partition_1_1_1_1}\n")
    
    # Final Summation
    total_count = num_partition_4 + num_partition_3_1 + num_partition_2_2 + num_partition_2_1_1 + num_partition_1_1_1_1
    
    print("--- Final Calculation ---")
    print("The total number is the sum of the counts for each partition:")
    print(f"Total = {num_partition_4} + {num_partition_3_1} + {num_partition_2_2} + {num_partition_2_1_1} + {num_partition_1_1_1_1} = {total_count}")

    print(f"\n<<< {total_count} >>>")

if __name__ == "__main__":
    main()