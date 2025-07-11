# The problem asks for the minimum integer n satisfying certain conditions
# based on the memory capacity of an agent.
#
# Let n be the length of the observation sequence.
# Let m be the number of states in the agent's memory (FSM).
#
# The agent can achieve a higher reward than a memoryless agent (m=1) if and only if
# it can configure its m-state FSM to distinguish between the observation sequence
# from corridor 1 (Omega_1) and corridor 2 (Omega_2).
#
# The problem asks for the minimum n such that there exists a pair of sequences
# (Omega_1, Omega_2) and a reward function for which:
# 1. An agent with m=2 states cannot get a higher reward than m=1.
#    This means NO 2-state FSM can distinguish Omega_1 and Omega_2.
# 2. An agent with m=3 states CAN get a higher reward than m=1.
#    This means there EXISTS a 3-state FSM that can distinguish Omega_1 and Omega_2.
#
# This is a known problem in automata theory: finding the length of the shortest
# pair of words that are k-indistinguishable but (k+1)-distinguishable. For k=2,
# we are looking for the minimal n for which there exist words that are
# 2-indistinguishable but 3-distinguishable.
#
# - For n=1, any two different sequences (e.g., (0) and (1)) are distinguishable by a 2-state FSM.
# - For n=2, it can be shown that for any pair of distinct length-2 sequences, a distinguishing 2-state FSM can be constructed.
# - For n=3, it is a known result in the field that there exist pairs of words that are indistinguishable by any 2-state FSM but are distinguishable by a 3-state FSM.
#
# Therefore, the minimum such length n is 3.

minimum_n = 3

print(minimum_n)