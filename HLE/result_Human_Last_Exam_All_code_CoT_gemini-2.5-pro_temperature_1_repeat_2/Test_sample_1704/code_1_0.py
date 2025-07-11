import math

def find_t4():
    """
    Finds the value of T(4), where T(n) is the sum of all elements in S(n),
    and S(n) is the set of all numbers that can be expressed as a sum of n
    distinct positive integers whose reciprocals sum to exactly 1.
    """
    n = 4
    sums = []
    
    # We are looking for x1, x2, x3, x4 such that:
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # with 1 < x1 < x2 < x3 < x4
    
    # From 1 < 4/x1, we get x1 < 4. So x1 can be 2 or 3.
    for x1 in range(2, 4):
        # From 1/x1 + 3/x2 > 1, we get x2 < 3*x1 / (x1 - 1)
        x2_upper_bound = math.ceil(3 * x1 / (x1 - 1))
        for x2 in range(x1 + 1, x2_upper_bound):
            # Let rem1 = 1 - 1/x1 - 1/x2.
            # From 1/x1 + 1/x2 + 2/x3 > 1, we get x3 < 2 / rem1.
            # To avoid float issues, we work with numerators and denominators.
            rem1_num = x1 * x2 - x2 - x1
            rem1_den = x1 * x2
            if rem1_num <= 0:
                continue
            
            x3_upper_bound = math.ceil((2 * rem1_den) / rem1_num)
            for x3 in range(x2 + 1, x3_upper_bound):
                # Let rem2 = 1 - 1/x1 - 1/x2 - 1/x3.
                # x4 = 1 / rem2.
                rem2_num = rem1_num * x3 - rem1_den
                rem2_den = rem1_den * x3
                
                if rem2_num <= 0:
                    continue
                
                # Check if x4 is an integer and if x4 > x3
                if rem2_den % rem2_num == 0:
                    x4 = rem2_den // rem2_num
                    if x4 > x3:
                        current_sum = x1 + x2 + x3 + x4
                        sums.append(current_sum)

    # Sort the sums to have a canonical representation for S(4)
    s4 = sorted(sums)
    total_sum_t4 = sum(s4)

    print(f"The elements of S(4) are the sums of each valid set of integers: {s4}")
    print("\nT(4) is the sum of these elements:")
    
    # Format the output string for the sum equation
    equation_str = " + ".join(map(str, s4))
    
    print(f"T(4) = {equation_str} = {total_sum_t4}")

find_t4()
<<<208>>>