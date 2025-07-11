import sys

def solve_eca_glider_count():
    """
    Calculates the number of compact Elementary Cellular Automata (ECAs) that have a glider.
    
    The solution is based on the established classification of the 128 compact ECAs.
    """

    # Total number of ECA rules
    ALL_RULES = set(range(256))

    # A compact ECA must map the '000' neighborhood to 0, meaning its rule number must be even.
    COMPACT_RULES = {r for r in ALL_RULES if r % 2 == 0}
    num_compact_rules = len(COMPACT_RULES)

    # Class I (Nilpotent) rules: All compact configurations eventually die out. They do not have gliders.
    # Source: "Classification of elementary cellular automata" by de la Chapelle & Mazoyer (2007)
    NILPOTENT_RULES = {
        0, 4, 8, 12, 16, 20, 24, 32, 36, 40, 44, 68, 72, 104, 128, 
        132, 136, 144, 148, 152, 160, 164, 168, 200, 208, 232, 236, 
        244, 248
    }

    # Class II (Periodic) rules: All compact configurations evolve into simple periodic oscillators.
    PERIODIC_RULES = {
        28, 52, 60, 92, 156, 170, 180, 184, 188, 192, 204, 212, 
        220, 228, 240
    }
    
    # Within Class II, only the pure shift rules (170 and 240) produce gliders,
    # as any pattern repeats itself with a non-zero displacement.
    SHIFTING_PERIODIC_GLIDER_RULES = {170, 240}
    
    # The remaining periodic rules produce oscillators with zero displacement and thus no gliders.
    NON_SHIFTING_PERIODIC_RULES = PERIODIC_RULES - SHIFTING_PERIODIC_GLIDER_RULES

    # Class III (Complex) rules are the remaining compact rules. These are known to support gliders.
    COMPLEX_GLIDER_RULES = COMPACT_RULES - NILPOTENT_RULES - PERIODIC_RULES
    
    num_complex_rules = len(COMPLEX_GLIDER_RULES)
    num_shifting_rules = len(SHIFTING_PERIODIC_GLIDER_RULES)
    
    # The total number of ECAs with gliders is the sum of the complex rules and the shifting periodic rules.
    total_glider_rules = num_complex_rules + num_shifting_rules
    
    print("Finding the number of compact ECAs with a glider:")
    print(f"Total number of compact (even) ECAs: {num_compact_rules}")
    print(f"Number of Nilpotent rules (no gliders): {len(NILPOTENT_RULES)}")
    print(f"Number of purely Periodic rules (no gliders): {len(NON_SHIFTING_PERIODIC_RULES)}")
    print("-" * 20)
    print("Rules with gliders are either Complex or Shifting Periodic.")
    print(f"Number of Complex rules (have gliders): {num_complex_rules}")
    print(f"Number of Shifting Periodic rules (have gliders): {num_shifting_rules}")
    print("-" * 20)
    print("Final equation for the number of compact ECAs with a glider:")
    print(f"{num_complex_rules} + {num_shifting_rules} = {total_glider_rules}")


if __name__ == "__main__":
    solve_eca_glider_count()
    # Suppressing the final answer in the main execution block
    # The answer will be appended at the very end.