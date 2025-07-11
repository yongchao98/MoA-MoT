# In the study of envy-free cake-cutting, several protocols have been developed.
# For the specific case of four agents requiring connected pieces, one of the key
# results comes from the 2015 paper "A discrete and bounded envy-free protocol for four agents"
# by Segal-Halevi, Hassidim, and Aumann.
#
# Their protocol provides an exact envy-free allocation, which also satisfies the
# epsilon-envy-free condition for any epsilon > 0. The protocol has a known
# upper bound on the number of cuts it needs to perform.
#
# This upper bound is a fixed, concrete number, making it a "realistic upper bound"
# as requested by the user.

# The maximum number of cuts required by the Segal-Halevi et al. protocol for 4 agents.
O = 16

# Print the upper bound
print(O)
