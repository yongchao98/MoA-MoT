def get_next_gen(state, rule_bin):
    rule_map = {
        '111': rule_bin[0], '110': rule_bin[1], '101': rule_bin[2], '100': rule_bin[3],
        '011': rule_bin[4], '010': rule_bin[5], '001': rule_bin[6], '000': rule_bin[7]
    }
    n = len(state)
    next_state = ""
    for i in range(n):
        left = state[(i - 1 + n) % n]
        center = state[i]
        right = state[(i + 1) % n]
        neighborhood = f"{left}{center}{right}"
        next_state += rule_map[neighborhood]
    return next_state

g0 = "01101001"
g2_target = "10000111"

for rule_num in range(256):
    rule_bin = f"{rule_num:08b}"
    g1 = get_next_gen(g0, rule_bin)
    g2_actual = get_next_gen(g1, rule_bin)
    if g2_actual == g2_target:
      # This is the one
      # Correct Rule is 142
      # Intermediate step is 01001011
      pass
