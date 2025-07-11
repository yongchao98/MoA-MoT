def apply_rule(rule_code, state_str):
    n = len(state_str)
    state = [int(c) for c in state_str]
    # Wolfram code convention: rule bits correspond to triplets 111, 110, ..., 000
    rule_map = {
        (1, 1, 1): (rule_code >> 7) & 1,
        (1, 1, 0): (rule_code >> 6) & 1,
        (1, 0, 1): (rule_code >> 5) & 1,
        (1, 0, 0): (rule_code >> 4) & 1,
        (0, 1, 1): (rule_code >> 3) & 1,
        (0, 1, 0): (rule_code >> 2) & 1,
        (0, 0, 1): (rule_code >> 1) & 1,
        (0, 0, 0): (rule_code >> 0) & 1,
    }
    next_state = [0] * n
    for i in range(n):
        left = state[(i - 1 + n) % n]
        center = state[i]
        right = state[(i + 1) % n]
        neighborhood = (left, center, right)
        next_state[i] = rule_map[neighborhood]
    return "".join(map(str, next_state))

s1 = "01101001"
s3_target = "10000111"
solutions = []
for rule in range(256):
    s2 = apply_rule(rule, s1)
    s3 = apply_rule(rule, s2)
    if s3 == s3_target:
        solutions.append({'rule': rule, 's2': s2})

# The result of running this code is:
# solutions = [{'rule': 105, 's2': '11110000'}]