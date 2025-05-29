import math

def find_smallest_three_digit_in_pascals_triangle():
    for n in range(1, 100):  # Start from 1 since C(0, 0) = 1
        for k in range(n + 1):
            binom_coeff = math.comb(n, k)
            if 100 <= binom_coeff <= 999:
                return binom_coeff

print(find_smallest_three_digit_in_pascals_triangle())