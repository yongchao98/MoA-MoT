import textwrap

def analyze_tree_ring_data():
    """
    Analyzes the factors influencing the declining 13C ratio in tree rings
    from 1886-1990 AD.
    """
    question = "Which of the factors below was the predominant factor influencing the declining 13C ratio observed in the tree ring data over this period?"

    options = {
        'A': 'An increase in tree ring thickness as the tree matures',
        'B': 'Periods of drought',
        'C': 'Increased photosynthetic reserves of starch fueling tree growth',
        'D': 'Thinning earlywood tree ring proportion',
        'E': 'Changes in the SE Asia monsoon'
    }

    # Print the question and the core scientific concept
    print(textwrap.fill("Task: " + question, 80))
    print("-" * 80)
    print(textwrap.fill("Scientific Background: A declining 13C ratio in tree rings indicates that the trees were experiencing less water stress over time. This allows the leaf stomata to stay open wider, increasing the uptake of the lighter 12C isotope relative to 13C during photosynthesis.", 80))
    print("-" * 80)
    print("Evaluating the options:\n")

    # Evaluate each option
    # A
    reasoning_a = "An increase in tree ring thickness is not a standard observation; rings typically get thinner as a tree ages. Furthermore, ring thickness is a measure of overall growth, not a direct cause of a change in isotopic ratios."
    print(f"Option A: {options['A']}")
    print(textwrap.fill(f"   - Analysis: {reasoning_a} This option is physiologically unlikely to be the predominant cause.\n", 80))

    # B
    reasoning_b = "Periods of drought cause water stress, which makes the tree close its stomata. This would lead to an INCREASE, not a decrease, in the 13C ratio because the tree is less able to discriminate against the heavier 13C isotope."
    print(f"Option B: {options['B']}")
    print(textwrap.fill(f"   - Analysis: {reasoning_b} This option would cause the opposite of the observed effect.\n", 80))

    # C
    reasoning_c = "While the use of stored energy (starch) can influence ring chemistry, it is not considered a primary driver of a consistent, century-long trend in isotopic composition compared to major atmospheric or climatic changes."
    print(f"Option C: {options['C']}")
    print(textwrap.fill(f"   - Analysis: {reasoning_c} This is a secondary effect, not a predominant factor.\n", 80))
    
    # D
    reasoning_d = "Changes in the proportion of earlywood (formed in spring) to latewood (formed in summer) are a result of changing growth conditions (like a shorter growing season or mid-summer drought), not the ultimate cause of the trend itself."
    print(f"Option D: {options['D']}")
    print(textwrap.fill(f"   - Analysis: {reasoning_d} This is a symptom of environmental change, not the root cause.\n", 80))

    # E
    reasoning_e = "The SE Asia monsoon is the primary source of precipitation for Chinese pine trees in many regions. A long-term trend towards a stronger or wetter monsoon would mean more water availability and higher humidity. This reduces water stress on the trees, leading to the observed decline in the 13C ratio."
    print(f"Option E: {options['E']}")
    print(textwrap.fill(f"   - Analysis: {reasoning_e} This is a plausible, large-scale climatic mechanism that directly explains the observed trend.\n", 80))

    # Conclusion
    print("-" * 80)
    print("Conclusion: Changes in a major climatic system like the SE Asia monsoon provide the most direct and powerful explanation for a long-term trend in regional water availability and the corresponding decline in the tree ring 13C ratio.")
    print("-" * 80)

# Execute the analysis
analyze_tree_ring_data()