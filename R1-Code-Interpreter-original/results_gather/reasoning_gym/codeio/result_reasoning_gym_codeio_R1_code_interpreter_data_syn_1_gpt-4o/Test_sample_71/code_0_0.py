def find_max_n_for_digit_sum(target_digit_sum):
    pre_numerator = 1
    cur_numerator = 2
    max_n = 2

    while True:
        temp = pre_numerator
        e_cont = 2 * max_n // 3 if max_n % 3 == 0 else 1
        pre_numerator = cur_numerator
        cur_numerator = e_cont * pre_numerator + temp

        digit_sum = sum(int(digit) for digit in str(cur_numerator))
        if digit_sum == target_digit_sum:
            return max_n

        max_n += 1

target_digit_sum = 346
max_n = find_max_n_for_digit_sum(target_digit_sum)
print(max_n)