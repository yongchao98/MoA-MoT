def get_bit(num, bit):
    return 1 if num & (1 << (bit - 1)) else 0

def set_bit(num, bit):
    return num | (1 << (bit - 1))

def clear_bit(num, bit):
    return num & ~(1 << (bit - 1))

def update_bit(num, bit, value):
    if value:
        return set_bit(num, bit)
    else:
        return clear_bit(num, bit)

def insertion(N, M, i, j):
    for bit in range(1, j - i + 2):
        b = get_bit(M, bit)
        N = update_bit(N, i + bit - 1, b)
    return N

def find_inputs_for_output(target_output):
    # Assume M is a small number, e.g., less than 256
    for M in range(0, 256):
        # Assume i and j are close to each other
        for i in range(1, 29):  # 28 is chosen to ensure j can be at most 32
            for j in range(i, min(i + 4, 33)):  # j is at most 4 bits away from i
                for N in range(0, 2**32):
                    if insertion(N, M, i, j) == target_output:
                        return N, M, i, j
    return None

# Target output
target_output = 1881421453

# Find inputs
result = find_inputs_for_output(target_output)
print(result)