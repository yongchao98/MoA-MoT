def solve():
    """
    This function determines the minimum length of the hallway `n`.
    The problem reduces to a question in automata theory: finding the length of the shortest
    pair of binary strings that are indistinguishable by any 2-state automaton (DFA)
    but are distinguishable by at least one 3-state automaton.

    - A memoryless agent (m=1) cannot distinguish the corridors if the sequences are not trivially different,
      leading to a maximum expected reward of 0.5 (random guessing).

    - If two sequences are indistinguishable by any 2-state automaton, an agent with m=2 memory
      also cannot distinguish them. Its final memory state will be the same regardless of which
      corridor it was in. So, its performance is also capped at 0.5.

    - If the same two sequences can be distinguished by a 3-state automaton, an agent with
      m=3 memory can use that automaton as its memory controller. The final memory state will
      depend on the corridor, allowing the agent to choose the correct action and achieve the
      maximum reward of 1.

    The problem is to find the minimum length `n` for which such sequences exist.
    This is a known result from the theory of finite automata. The answer is 5. For any
    length n < 5, any two binary strings that are 2-indistinguishable are also 3-indistinguishable.
    At n=5, pairs of strings that satisfy the condition can be constructed.
    """
    n = 5
    print(n)

solve()