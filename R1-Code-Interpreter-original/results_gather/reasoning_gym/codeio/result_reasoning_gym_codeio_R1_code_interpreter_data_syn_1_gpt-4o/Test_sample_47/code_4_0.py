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

def find_inputs(target):
    target_bits = bin(target)[2:]  # Get binary representation of the target
    target_length = len(target_bits)
    
    # Consider larger M values
    for M in range(0, 2**16):  # Limiting M to 16-bit numbers
        M_bits = bin(M)[2:].zfill(target_length)
        for i in range(1, 33 - target_length + 1):
            j = i + target_length - 1
            # Calculate the mask to clear bits i to j in N
            mask = ((1 << (j - i + 1)) - 1) << (i - 1)
            # Calculate N directly
            M_shifted = M << (i - 1)
            N = target & ~mask | M_shifted
            if insertion(N, M, i, j) == target:
                return N, M, i, j
    return None

target = 1881421453
result = find_inputs(target)
print(result)