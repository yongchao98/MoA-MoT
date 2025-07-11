import numpy as np

# This is a simulator for the provided RASPy-like code.
# We implement the semantics of the operations as Python functions.

# Helper to create numpy arrays from sequences
def to_arr(seq):
    if isinstance(seq, list) and len(seq) > 0 and isinstance(seq[0], str):
        return np.array(list(seq))
    if isinstance(seq, str):
        return np.array(list(seq))
    return np.array(seq)

# Simulating RASPy operations
def raspy_where(cond, a, b):
    cond, a, b = to_arr(cond), to_arr(a), to_arr(b)
    # Ensure a and b are broadcastable to the condition's shape
    if a.shape != cond.shape: a = np.broadcast_to(a, cond.shape)
    if b.shape != cond.shape: b = np.broadcast_to(b, cond.shape)
    return np.where(cond, a, b)

def raspy_map(seq, func):
    return np.array([func(x) for x in to_arr(seq)])

def raspy_d(seq):
    # Cumulative sum
    return np.cumsum(to_arr(seq))

def raspy_f(i, default, seq):
    seq = to_arr(seq)
    res = np.full_like(seq, fill_value=default)
    if i > 0:
        res[i:] = seq[:-i]
    elif i < 0:
        res[:i] = seq[-i:]
    else: # i == 0
        res = seq
    return res

# --- Implementing the provided functions ---

def a(seq):
    return raspy_map(seq, lambda x: ord(x) - ord('0'))

def l(default, sop):
    sop_arr = to_arr(sop)
    c = np.sum(sop_arr == '_')
    # This is a right shift
    return raspy_f(c, default, sop_arr)
    
def m(v, i, sop, default="0"):
    sop_arr = to_arr(sop)
    split_indices = np.where(sop_arr == v)[0]
    split_point = split_indices[0] if len(split_indices) > 0 else -1
    
    indices = np.arange(len(sop_arr))
    if i: # True, part before split
        seq_to_l = raspy_where(indices < split_point, sop_arr, "_")
        return l(default, seq_to_l)
    else: # False, part after split
        return raspy_where(indices > split_point, sop_arr, default)

def n(match, seq):
    match, seq = to_arr(match), to_arr(seq)
    x = raspy_d(raspy_map(match, int))
    y = np.full_like(seq, fill_value=seq[-1]) # Default fill for y
    
    # Right-to-left propagation
    for i in range(len(seq) - 2, -1, -1):
        if not match[i]:
            # Find first non-match to the right
            k = i + 1
            while k < len(seq) and not match[k]:
                k +=1
            if k < len(seq):
                y[i] = seq[k]
            else: # No non-match found to the right
                y[i] = seq[i] # Or a default, behavior is ambiguous, but this seems safe
        else:
            y[i] = seq[i]

    # This simplified loop implements the described right-to-left propagation logic
    # which is what the complex RASPy `y` calculation effectively does.
    result = np.copy(seq)
    for i in range(len(seq) - 2, -1, -1):
        if result[i] == '<' and result[i+1] == '1':
            result[i] = '1'
    result = raspy_where(result == '<', '0', result)
    return result

def calculate_sum(sop):
    # This simulates the logic of aa, bb, cc in function v
    sop = str(sop)
    # 1. Parse and Align
    num1_aligned_str = m("+", True, sop, default="0")
    num2_aligned_str = m("+", False, sop, default="0")

    num1_int = a(num1_aligned_str)
    num2_int = a(num2_aligned_str)
    
    # 2. Digit-wise Sum
    aa = num1_int + num2_int
    
    # 3. Carry Calculation
    carry_map = raspy_where(aa > 9, "1", raspy_where(aa == 9, "<", "0"))
    
    # Left shift carry
    bb = raspy_f(-1, "0", carry_map)

    # 4. Carry Propagation
    carry_final_str = n(bb != '<', bb)
    carry_final = a(carry_final_str)

    # 5. Final Sum
    cc = (aa + carry_final) % 10
    
    # 6. Formatting
    result_str = "".join(map(str, cc))
    # Omitting leading zeros
    return result_str.lstrip('0') or '0'


def solve():
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = calculate_sum(input1)
    output2 = calculate_sum(input2)

    print(f"{output1};{output2}")

solve()