def get_second_distinct(chunk):
    """Finds the second distinct element in a list."""
    seen = []
    for x in chunk:
        if x not in seen:
            seen.append(x)
        if len(seen) == 2:
            return seen[1]
    return None

def generate_sequence():
    """Generates the sequence based on the chunk rules."""
    
    # History of chunk properties
    s_history = []
    len_history = []
    full_sequence = []
    
    # --- k=1 ---
    k = 1
    s_k = 3
    d_k = 'down'
    s_history.append(s_k)
    
    # Down rule for k=1
    c_k = list(range(s_k, 0, -1))
    len_history.append(len(c_k))
    full_sequence.extend(c_k)
    
    # --- k=2 ---
    k = 2
    s_k = get_second_distinct(c_k)
    d_k = 'up'
    s_history.append(s_k)
    s_prev = s_history[-2]
    
    # Up rule for k<=2
    c_k = list(range(s_k, s_prev + 1))
    len_history.append(len(c_k))
    full_sequence.extend(c_k)

    # --- k=3 ---
    k = 3
    s_k = get_second_distinct(c_k)
    d_k = 'down'
    s_history.append(s_k)
    len_prev = len_history[-1]
    
    # Down rule for k>1
    c_k = [s_k] * len_prev + [s_k - 1] * len_prev
    len_history.append(len(c_k))
    full_sequence.extend(c_k)
    
    # --- k=4 (Prediction) ---
    k = 4
    s_k = get_second_distinct(c_k)
    d_k = 'up'
    s_history.append(s_k)
    s_prev = s_history[-2]

    # Up rule for k>2
    part1 = list(range(s_k, s_prev + 1))
    # Count down from previous S, but exclude the first element
    part2 = list(range(s_prev - 1, 0, -1))
    c_k_next = part1 + part2
    
    print(f"The given sequence is: {' '.join(map(str, full_sequence))}")
    # The final equation part is asking to output each number clearly.
    print(f"The next 4 elements are: {c_k_next[0]} {c_k_next[1]} {c_k_next[2]} {c_k_next[3]}")

generate_sequence()