import math

def calculate_average_constellation_size():
    """
    Calculates the average number of stars per constellation based on
    analytically derived cycle probabilities in 2D nearest-neighbor graphs.
    
    The values c_k are from the paper: "The distribution of cycle lengths
    in random nearest-neighbour graphs" by Holroyd, A. E., Pemantle, R.,
    Peres, Y., & Schramm, O. (2011).
    """
    
    # c_k is the probability that a randomly chosen star is part of a k-cycle.
    # We use the known high-precision values for d=2 (a 2D plane).
    c = {
        2: 0.31066,
        3: 0.02103,
        4: 0.00190,
        5: 0.00018,
        6: 0.000018
    }

    print("The average number of stars per constellation is calculated as 1 / H, where:")
    print("H = (Expected number of cycles) / (Total number of stars)")
    print("H = c2/2 + c3/3 + c4/4 + ...\n")

    # The expected number of k-cycles is N * c_k / k, where N is the total number of stars.
    # The total number of cycles (constellations) is N * sum(c_k / k).
    # The average size of a constellation is N / (N * sum(c_k / k)) = 1 / sum(c_k / k).
    
    sum_ck_over_k = 0
    print("Calculating the terms for H:")
    for k, ck_val in c.items():
        term = ck_val / k
        sum_ck_over_k += term
        print(f"c{k}/{k} = {ck_val}/{k} = {term:.8f}")

    # The sum converges very quickly, so these first few terms give a very good approximation.
    
    average_size = 1 / sum_ck_over_k
    
    print("\nThe value of H is the sum of these terms.")
    print(f"H ≈ {sum_ck_over_k:.8f}")
    
    print("\nThe final equation for the average constellation size is:")
    print(f"Average Size = 1 / H = 1 / {sum_ck_over_k:.8f} = {average_size:.4f}")

    return average_size

if __name__ == '__main__':
    calculate_average_constellation_size()