import functools
import operator

def solve_nim_variant():
    """
    Solves the Nim variant game for five scenarios and determines the winner.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    def h(a):
        """
        Calculates the Grundy value for a pile of size 'a'
        when the Zero-Move has not been used.
        """
        if a == 0:
            return 0
        if a % 2 == 0:  # a is even and positive
            return a - 1
        else:  # a is odd
            return a + 1

    result_string = ""
    for i, piles in enumerate(scenarios):
        print(f"Scenario ({i+1}) with piles {piles}:")
        
        grundy_values = [h(p) for p in piles]
        print(f"Grundy values: {grundy_values}")

        nim_sum = functools.reduce(operator.xor, grundy_values)
        
        equation_str = " ^ ".join(map(str, grundy_values))
        print(f"Nim-sum: {equation_str} = {nim_sum}")
        
        if nim_sum == 0:
            winner = "Bob"
            result_string += "B"
        else:
            winner = "Alice"
            result_string += "A"
            
        print(f"Winner: {winner}")
        print("-" * 20)

    print(f"Final answer string: {result_string}")
    print(f"<<<{result_string}>>>")

solve_nim_variant()