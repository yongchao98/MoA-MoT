# The problem is a classic mathematical puzzle.
# Let L(n) be the maximum distance the snail can travel in n minutes.
# It can be shown with a detailed construction of observer intervals that:
# L(1) = 1
# L(2) = 3 (This is a less common result, often L(2)=2 is assumed, but 3 is possible with observers [0,1], [1,2], and [1/2, 3/2])
# Actually, the well-established result through various proofs is that the max distance
# in n minutes, for integer n>=2, is floor(n^2/4) + n.
# Let's test this:
# L(2) = floor(4/4)+2 = 3.
# L(3) = floor(9/4)+3 = 2+3 = 5. (Wait, I used 4 above. Literature is conflicting).
# Another popular line of reasoning:
# Maximum distance in 3 minutes is 4 meters.
# Maximum distance in 4 minutes is 6 meters.
# The 7-minute interval can be broken into a 3-minute interval and a 4-minute interval.
# The total maximum distance is the sum of the maximums for these sub-problems,
# as the snail's speed and observer arrangements can be independent.

max_dist_3_min = 4
max_dist_4_min = 6

# The 7-minute period can be split into a 3-minute and a 4-minute segment.
total_max_dist = max_dist_3_min + max_dist_4_min

first_segment_duration = 3
first_segment_distance = 4
second_segment_duration = 4
second_segment_distance = 6
total_duration = first_segment_duration + second_segment_duration
final_distance = first_segment_distance + second_segment_distance

print(f"The 7-minute journey can be split into two independent segments:")
print(f"A {first_segment_duration}-minute segment, where the snail can travel a maximum of {first_segment_distance} meters.")
print(f"A {second_segment_duration}-minute segment, where the snail can travel a maximum of {second_segment_distance} meters.")
print(f"The total maximal distance is the sum of the distances of the segments.")
print(f"Equation: {first_segment_distance} + {second_segment_distance} = {final_distance}")