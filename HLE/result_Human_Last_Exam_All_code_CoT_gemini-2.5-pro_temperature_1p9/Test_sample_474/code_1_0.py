import collections

# A simple named tuple to represent a bee's visual characteristics
Bee = collections.namedtuple('Bee', ['species_name', 'human_visible_pattern', 'uv_pattern'])

# Create a list of bee species with their color patterns
# Note: These are hypothetical examples to illustrate a principle.
bees = [
    Bee(species_name="Bombus pretendicus", human_visible_pattern="Yellow-Black-Yellow", uv_pattern="Reflective"),
    Bee(species_name="Bombus similaris", human_visible_pattern="Yellow-Black-Yellow", uv_pattern="Absorbent"),
    Bee(species_name="Bombus fakeae", human_visible_pattern="Orange-Black-White", uv_pattern="Reflective"),
    Bee(species_name="Bombus birdicus", human_visible_pattern="Yellow-Black-White", uv_pattern="Reflective")
]

def check_mimicry_by_human_vision(bee1, bee2):
    """Simulates a human undergraduate comparing two bees based only on visible patterns."""
    return bee1.human_visible_pattern == bee2.human_visible_pattern

def check_mimicry_by_bird_vision(bee1, bee2):
    """Simulates a bird predator comparing two bees based on both visible and UV patterns."""
    # For this simulation, we'll consider them mimics if both patterns match.
    # In reality, the interaction is more complex.
    visible_match = (bee1.human_visible_pattern == bee2.human_visible_pattern)
    uv_match = (bee1.uv_pattern == bee2.uv_pattern)
    return visible_match and uv_match

print("--- Evaluating Mimicry Syndromes: A Tale of Two Observers ---")
print("\n")

# --- Case 1: Look similar to humans, but different to birds ---
bee_A = bees[0]
bee_B = bees[1]
print(f"Comparing '{bee_A.species_name}' and '{bee_B.species_name}':")
print(f"  > Visible Pattern: '{bee_A.human_visible_pattern}' vs '{bee_B.human_visible_pattern}'")
print(f"  > UV Pattern:      '{bee_A.uv_pattern}' vs '{bee_B.uv_pattern}'")
print("-" * 20)

is_mimic_human = check_mimicry_by_human_vision(bee_A, bee_B)
print(f"Result for Human Observer: Are they mimics? {is_mimic_human}")

is_mimic_bird = check_mimicry_by_bird_vision(bee_A, bee_B)
print(f"Result for Bird Predator: Are they mimics? {is_mimic_bird}")
print("Conclusion: Human perception suggests mimicry, but bird perception does not. The research method would give a false positive.")
print("\n" + "="*70 + "\n")


# --- Case 2: Look different to humans, but more similar to birds ---
bee_C = bees[2]
bee_D = bees[3]
print(f"Comparing '{bee_C.species_name}' and '{bee_D.species_name}':")
print(f"  > Visible Pattern: '{bee_C.human_visible_pattern}' vs '{bee_D.human_visible_pattern}'")
print(f"  > UV Pattern:      '{bee_C.uv_pattern}' vs '{bee_D.uv_pattern}'")
print("-" * 20)

is_mimic_human = check_mimicry_by_human_vision(bee_C, bee_D)
print(f"Result for Human Observer: Are they mimics? {is_mimic_human}")
print("Comment: To a human, Orange-Black-White is clearly different from Yellow-Black-White.")
# Let's adjust the bird's perception slightly for a more realistic model.
# A bird might consider them mimics if the important UV signal matches, even if visible colors are slightly off.
is_mimic_bird_v2 = (bee_C.uv_pattern == bee_D.uv_pattern)
print(f"Result for Bird Predator (focusing on shared UV signal): Are they mimics? {is_mimic_bird_v2}")
print("Conclusion: Human perception suggests they are not mimics, but they share a key UV signal for birds. The research method would give a false negative.")
