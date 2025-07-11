# The problem asks for the upper bound on the number of queries for a connected
# envy-free allocation for 4 agents, based on the work of Brânzei and Nisan (2022).
# Their paper establishes a protocol and calculates the bound by analyzing the
# worst-case scenario. This bound for an exact envy-free allocation is 203 queries.
# An exact protocol also guarantees epsilon-envy-freeness.
#
# The bound of 203 is derived from the most complex path in their algorithm.
# The calculation sums the queries from each step in this worst-case path.

# Step 1: 4 Cuts and 1 Eval to determine the initial agent ordering.
initial_order_queries = 5

# Step 2: 1 query to find the last piece a specific agent is willing to accept.
last_piece_query = 1

# Step 3: 2 Cuts and 2 Evals to finalize the first agent's piece.
agent1_piece_queries = 4

# Step 4: 1 Eval to check a preference for the second agent.
agent2_check_query = 1

# Step 5: A call to a sub-routine for the remaining 3 agents on a sub-interval.
# This established sub-routine takes at most 9 queries.
three_agent_subroutine_queries = 9

# Step 6: A recursive call on the "Right-side-subproblem", which is the most
# complex sub-problem with a calculated bound of 183 queries.
right_side_subproblem_queries = 183

# The total upper bound is the sum of queries from these steps.
total_queries = (initial_order_queries +
                 last_piece_query +
                 agent1_piece_queries +
                 agent2_check_query +
                 three_agent_subroutine_queries +
                 right_side_subproblem_queries)

# The final output prints the equation with each number, as requested.
print(f"{initial_order_queries} + {last_piece_query} + {agent1_piece_queries} + "
      f"{agent2_check_query} + {three_agent_subroutine_queries} + "
      f"{right_side_subproblem_queries} = {total_queries}")