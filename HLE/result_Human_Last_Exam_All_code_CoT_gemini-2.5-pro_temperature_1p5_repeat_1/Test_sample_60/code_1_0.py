import collections

def parse_rules(rule_string):
    """Parses a string representation of TM rules into a dictionary."""
    parts = rule_string.split()
    states = 'ABCDE'
    rules = {}
    if len(parts) != len(states) * 2:
        raise ValueError("Invalid number of rules in string")
    
    for i, state_char in enumerate(states):
        rule_0_str = parts[i * 2]
        rule_1_str = parts[i * 2 + 1]
        rules[state_char] = {
            0: {'write': int(rule_0_str[1]), 'move': -1 if rule_0_str[2] == 'L' else 1, 'next_state': rule_0_str[0]},
            1: {'write': int(rule_1_str[1]), 'move': -1 if rule_1_str[2] == 'L' else 1, 'next_state': rule_1_str[0]}
        }
    return rules

def simulate_tm(rules, max_steps=100000):
    """Simulates a Turing Machine and returns the number of steps to halt."""
    tape = collections.defaultdict(int)
    head_pos = 0
    state = 'A'
    steps = 0
    while state != 'H':
        if steps >= max_steps:
            return -1  # Indicates it did not halt within the limit
        current_symbol = tape[head_pos]
        rule = rules[state][current_symbol]
        tape[head_pos] = rule['write']
        head_pos += rule['move']
        state = rule['next_state']
        steps += 1
    return steps

def solve():
    """
    Simulates three Turing Machines and identifies the one that halts after the most steps.
    """
    machines = {
        1: "B1L B0R C1R E0L A1L D0R C0R A1R C1L H1L",
        2: "B1L C1R C1L B1L D1L E0R A1R D1R H1L A0R",
        3: "B1L A1L C1R B1R A1L D1R A1L E1R H1L C0R"
    }
    
    results = {}
    for i, machine_str in machines.items():
        rules = parse_rules(machine_str)
        steps = simulate_tm(rules)
        results[i] = steps

    winner_id = -1
    max_steps = -1
    for i, steps in results.items():
        print(f"Machine {i} halted after {steps} steps.")
        if steps > max_steps:
            max_steps = steps
            winner_id = i

    print(f"\nThe machine that halts after the most number of steps is Machine {winner_id} with {max_steps} steps.")

if __name__ == '__main__':
    solve()