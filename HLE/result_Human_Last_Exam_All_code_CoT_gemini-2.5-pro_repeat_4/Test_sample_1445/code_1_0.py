# Plan:
# The problem asks for the maximum of the minimum operations required to transform any sequence S to any sequence T.
# Let's consider the transformation between a uniform sequence and a maximally alternating sequence, as this represents a "maximal change" in structure.
#
# Let S be a uniform sequence of 100 '0's.
# Let T be a maximally alternating sequence starting with '1', i.e., '1010...10'.
#
# We can calculate the number of operations by considering an intermediate path through a uniform sequence.
#
# Path 1: S -> '00...0' -> T
# Cost(S -> '00...0'): S is already '00...0', so the cost is 0.
# Cost('00...0' -> T):
#   - T has 50 '1'-blocks and 50 '0'-blocks.
#   - To build T from '00...0', we must insert the 50 '1'-blocks.
#   - Inserting 50 '1'-blocks into '00...0' results in the sequence '0101...01' (let's call it T_prime). Cost = 50 ops.
#   - Now we must transform T_prime ('0101...01') into T ('1010...10'). This requires removing the first '0'-block and adding a '0'-block at the end. Cost = 2 ops.
#   - Total cost for '00...0' -> T is 50 + 2 = 52.
# Total cost for Path 1 = 0 + 52 = 52.
#
# Path 2: S -> '11...1' -> T
# Cost(S -> '11...1'):
#   - S is '00...0'. We remove the '0'-block (1 op) and insert a '1'-block (1 op). Cost = 2.
# Cost('11...1' -> T):
#   - T has 50 '0'-blocks. To build T from '11...1', we must insert these 50 '0'-blocks.
#   - This correctly constructs '1010...10'. Cost = 50 ops.
# Total cost for Path 2 = 2 + 50 = 52.
#
# Both paths give a cost of 52. This is the cost for this specific S and T pair and represents the worst-case scenario.

cost_path1_step1 = 0
cost_build_t_prime = 50
cost_transform_t_prime_to_t = 2
cost_path1_step2 = cost_build_t_prime + cost_transform_t_prime_to_t
total_cost_path1 = cost_path1_step1 + cost_path1_step2

cost_path2_step1 = 2
cost_path2_step2 = 50
total_cost_path2 = cost_path2_step1 + cost_path2_step2

# The minimum number of operations 'n' is the maximum of these costs over all possible S and T.
# The calculated value of 52 for this "most different" pair is the maximum.
n = max(total_cost_path1, total_cost_path2)

print(f"The cost of Path 1 is {cost_path1_step1} + ({cost_build_t_prime} + {cost_transform_t_prime_to_t}) = {total_cost_path1}")
print(f"The cost of Path 2 is {cost_path2_step1} + {cost_path2_step2} = {total_cost_path2}")
print(f"The minimum number of operations n needed is {n}")
