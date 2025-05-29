def is_happy(num):
    seen = set()
    while num != 1 and num not in seen:
        seen.add(num)
        num = sum(int(digit) ** 2 for digit in str(num))
    return num == 1

def find_start_num(target_happy_numbers):
    start_num = 1
    while True:
        happy_numbers = []
        num = start_num
        while len(happy_numbers) < 8:
            if is_happy(num):
                happy_numbers.append(num)
            num += 1
        if happy_numbers == target_happy_numbers:
            return start_num
        start_num += 1

target_happy_numbers = [464, 469, 478, 487, 490, 496, 536, 556]
start_num = find_start_num(target_happy_numbers)
print(start_num)