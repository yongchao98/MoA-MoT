import numpy as np

# This script simulates the provided RASPy code to solve the problem.
# The user's prompt has slightly conflicting formatting instructions.
# I will adhere to the specific "output1;output2" format requested for the final answer,
# as it is more precise than the general instruction about "the final equation".

# --- Python implementation of the RASPy functions ---

def as_seq(val, length):
    """Helper to handle scalar broadcasting like in RASPy."""
    if isinstance(val, (list, np.ndarray)):
        return np.array(val)
    return np.full(length, val)

def where(cond, true_val, false_val):
    """Numpy's equivalent of where."""
    return np.where(cond, true_val, false_val)

def a(seq):
    """Converts a sequence of digit characters to integers."""
    return np.array([ord(x) - ord('0') if x.isdigit() else 0 for x in seq], dtype=int)

def c(seq):
    """Converts a sequence of ASCII codes to a string."""
    return "".join([chr(int(x)) for x in seq])

def d(seq):
    """Calculates the cumulative sum of a sequence."""
    return np.cumsum(as_seq(seq, len(seq)))

def f(i, default, seq):
    """Shifts a sequence right by i positions, or left if i is negative."""
    seq = np.array(list(seq))
    res = np.full_like(seq, fill_value=default, dtype=seq.dtype)
    if i > 0:  # right shift
        if len(seq) > i:
            res[i:] = seq[:-i]
    elif i < 0:  # left shift
        i = -i
        if len(seq) > i:
            res[:-i] = seq[i:]
    else:  # i == 0
        res = seq
    return res

def h(i, default, seq):
    """Left shift by i-1."""
    return f(-(i - 1), default, seq)

def i_func(i_val, default, seq):
    """A combination of a right shift and a left shift."""
    seq = np.array(list(seq))
    # Right shift by (i_val - 3)
    seq_1 = f(i_val - 3, default, seq)
    # Left shift by (i_val - 3)
    seq_2 = f(-(i_val - 3), default, seq_1)
    return seq_2

def j(seq):
    """Finds the minimum value in a sequence and broadcasts it."""
    # Analysis showed this function effectively calculates the minimum value.
    min_val = np.min(seq[np.where(seq != '_')]) if '_' in seq else np.min(seq)
    return np.full(len(seq), min_val)

def l(default, sop):
    """Right-aligns content in a sequence by shifting right by the count of '_'."""
    sop = np.array(list(sop))
    c_val = np.sum(sop == '_')
    return f(c_val, default, sop)

def m(v_char, i_bool, sop, default="0"):
    """Splits the sop string at v_char and aligns the parts."""
    sop_list = list(sop)
    length = len(sop_list)
    indices = np.arange(length)
    try:
        split_point = sop_list.index(v_char)
    except ValueError:
        split_point = -1

    if i_bool:  # Left part of the '+'
        w1 = where(indices < split_point, sop_list, "_")
        return l(default, w1)
    else:  # Right part of the '+'
        return where(indices > split_point, sop_list, default)

def n(match, seq):
    """Propagates carry values ('0' or '1') across placeholder ('<') characters."""
    result = np.array(list(seq))
    for i in range(len(result) - 2, -1, -1):
        if result[i] == '<':
            result[i] = result[i+1]
    return result

def s(sop):
    """A complex function that was analyzed to be equivalent to calculating 3 * (number of '7's)."""
    sop = np.array(list(sop))
    n7 = np.sum(sop == '7')
    return np.full(len(sop), 3 * n7)

def t(seq):
    """Left-aligns a sequence by removing leading '_' characters through shifting."""
    seq = np.array(list(seq))
    try:
        first_char_idx = np.min(np.where(seq != '_'))
        return f(-first_char_idx, '_', seq)
    except ValueError:  # Handles case where sequence is all '_'
        return seq

def u(sop):
    """The trigger function. Returns [1,1,...] if the '777' condition is met, else [0,0,...]."""
    length = len(sop)
    indices = np.arange(length)

    # Split sop into 3 parts: first 3, middle 3, last 3
    part1__ = where(indices < 3, list(sop), '_')
    part2_ = where((indices > 2) & (indices < 6), list(sop), '_')
    part3_ = where(indices > length - 4, list(sop), '_')
    
    # Left-align them
    part1 = part1_ # Already aligned
    part2 = t(part2_)
    part3 = t(part3_)

    # Get the '7' count metric for each part
    s_p1 = s(part1)
    s_p2 = s(part2)
    s_p3 = s(part3)
    
    # Combine the metrics into a single sequence [3*N7_p1, 3*N7_p2, 3*N7_p3, 0, ...]
    h_p1 = h(length, 0, s_p1)
    h_p2 = h(length, 0, s_p2)
    k_p2 = f(1, 0, h_p2)
    h_p3 = h(length, 0, s_p3)
    n_p3 = f(2, 0, h_p3)
    oo = h_p1 + k_p2 + n_p3

    # Isolate the metrics and check if all are >= 1
    pp = i_func(length, 1, oo)
    qq = j(pp)
    return qq

def v(sop):
    """The main function."""
    length = len(sop)
    indices = np.arange(length)

    # --- Path 2: The trigger check ---
    dd = u(sop) # This will be [1,1,...] if triggered, [0,0,...] otherwise

    # If not triggered, calculate the sum
    if dd[0] == 0:
        num1_padded = m("+", True, sop, default="0")
        num2_padded = m("+", False, sop, default="0")
        
        aa = a(num1_padded) + a(num2_padded)
        carry_info = where(aa > 9, "1", where(aa == 9, "<", "0"))
        shifted_carry_info = f(-1, "0", carry_info)
        propagated_carries = n(shifted_carry_info != "<", shifted_carry_info)
        bb = a(propagated_carries)
        
        cc = (aa + bb) % 10
        # Convert digits to string, removing padding
        result_str = "".join(map(str, cc)).lstrip('0')
        return result_str if result_str else "0"

    # If triggered, compute the "pwned" message
    else:
        ee = np.array([103, 101, 116, 32, 112, 119, 110, 101, 100] + [33] * (length))
        ff = ee[:length]
        
        # Aesthetic change: last char becomes '1' if index > 10
        if length - 1 > 10:
            ff[length - 1] = 49  # ASCII for '1'
        
        return c(ff)

def solve():
    """Processes the inputs and prints the final result."""
    input1 = "734107+4295754"
    input2 = "5429141+142196"

    output1 = v(input1)
    output2 = v(input2)
    
    print(f"{output1};{output2}")

solve()