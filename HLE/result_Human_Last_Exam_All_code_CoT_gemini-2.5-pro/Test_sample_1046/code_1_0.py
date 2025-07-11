# The user wants to understand how a bacterium with a stable genome (Bacterium 2)
# can acquire drug resistance as quickly as a bacterium that frequently uses
# lateral gene transfer (Bacterium 1).

# Step 1: Define the problem.
# Bacterium 1 uses a high-frequency method (lateral transfer) to gain resistance.
# Bacterium 2 uses a low-frequency method (rare mutations) to gain resistance.
# The paradox is that they achieve resistance at an equal pace.

# Step 2: Analyze the requirements for Bacterium 2 to keep pace.
# For a low-frequency event to have the same overall effect as a high-frequency one,
# the impact of that single low-frequency event must be extraordinarily large.

# Step 3: Evaluate the options based on their impact.
# - A simple resistance mutation is a small-impact event.
# - A resistance mutation with a compensatory mutation is a medium-impact event (it makes the strain viable and competitive).
# - A resistance mutation that confers cross-resistance (resistance to multiple drugs) is a high-impact event.
# - A cross-resistance mutation that is ALSO accompanied by compensatory mutations is the highest-impact scenario.
#   This creates a "super strain" that is resistant to multiple drugs and has no fitness cost,
#   allowing it to proliferate extremely rapidly under selection pressure.

# Step 4: Conclude the best explanation.
# Option B describes this highest-impact scenario: the acquisition of compensatory mutations
# to increase fitness, combined with the powerful effect of cross-resistance. This provides
# the most complete explanation for how the second bacterial population could keep pace.

final_explanation = (
    "The correct answer is B. "
    "For the second bacterium, which relies on rare mutations, to match the pace of the first, "
    "a rare evolutionary event must have a massive payoff. "
    "Option B provides two key mechanisms that, together, create this massive payoff:\n"
    "1. Cross-resistance: A single mutation provides resistance to multiple drugs, a huge evolutionary shortcut.\n"
    "2. Compensatory mutations: These eliminate the fitness cost often associated with resistance, "
    "allowing the new multi-drug resistant strain to spread rapidly through the population.\n"
    "This combination creates a highly fit 'super strain' that can proliferate as quickly as resistance spreads "
    "via high-frequency lateral transfer in the first bacterium."
)

print(final_explanation)